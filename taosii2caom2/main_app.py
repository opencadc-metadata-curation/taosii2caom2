# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

"""
ML: 22-09-20 - teleconference:
HDF5 arrangement: four axes: epoch (time), telescope, x, y

What I did:
Time:

# mjd start

>>> t = Time(1561057509.599084, format='unix')
>>> t.format = 'mjd'
>>> t.value
58654.795249989395
>>> 58654.795249989395
58654.795249989395

# mjd end

>>> t = Time(0.05000000074505806+1561057509.599084, format='unix')
>>> t.format = 'mjd'
>>> t.value
58654.7952505681

start = RefCoord(0.5, mjd_start)
end = RefCoord(1.5, mjd_end)


telescope location:

>>> from astropy.coordinates import EarthLocation
>>> x = EarthLocation.of_site('spm')
Downloading http://data.astropy.org/coordinates/sites.json
|=============================================================================|  23k/ 23k (100.00%)         0s
>>> x
<EarthLocation (-2354953.99637757, -4940160.36363812, 3270123.70695983) m>
>>>


energy:

min = 400 nm
max = 800 nm
central = min + max / 2.0 = 1200.0 / 2.0 == 600.0 nm
fwhm = 400 nm
ref_coord1 = RefCoord(0.5, central_wl - fwhm / 2.0) == 100.0 nm
ref_coord2 = RefCoord(1.5, central_wl + fwhm / 2.0) == 500.0 nm

Position:

# assume equatorial coordinates

>>> from astropy import wcs
>>> w = wcs.WCS(naxis=2)
>>> w.wcs.crpix = [-3850, 2310]
>>> w.wcs.crval = [72.0, 20.77222]
>>> w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
>>> print(w.to_header())

"""

import h5py
import logging

from astropy import wcs
from datetime import datetime
from os.path import basename

from caom2 import CalibrationLevel, ProductType, DataProductType
from caom2 import ObservationIntentType

from caom2utils.caom2blueprint import Hdf5ObsBlueprint, Hdf5Parser
from caom2pipe.astro_composable import add_as_s, build_ra_dec_as_deg
from caom2pipe.caom_composable import Fits2caom2Visitor, TelescopeMapping
from caom2pipe.manage_composable import CadcException, StorageName, to_float, ValueRepairCache


__all__ = ['TAOSII2caom2Visitor', 'TAOSIIName']


class TaosiiValueRepair(ValueRepairCache):

    VALUE_REPAIR = {
        'chunk.position.axis.axis1.ctype': {
            'RA--TAN-SIP': 'RA---TAN-SIP',
            'DEC-TAN_SIP': 'DEC--TAN-SIP',
        },
    }

    def __init__(self):
        self._value_repair = TaosiiValueRepair.VALUE_REPAIR
        self._key = None
        self._values = None
        self._logger = logging.getLogger(self.__class__.__name__)


class BasicMapping(TelescopeMapping):

    def __init__(self, storage_name, h5file, clients, observable, observation):
        super().__init__(storage_name, headers=None, clients=clients, observable=observable, observation=observation)
        self._h5_file = h5file
        self._prefix = None

    def accumulate_blueprint(self, bp):
        """
        JJK - 17-03-22

         v = {"Observation":
             {"observationID":  f"{fobj['header'][]}",
              "type": f"{fobj['header']['run']['obstype']}",
              "metaRelease": now,
              "algorithm": {"name": "/header/run/exptype"},
              "instrument": {"name": "/header/run/origin",
                             "keywords": f"{fobj['header']['run']['exptype']}"},
              "proposaal": {"id": "/header/run/run_seq",
                            "pi": f"{fobj[['header']['run']['observer']}",
                            "project": "TAOS2",
                            "title": "Transneptunian Automated Occultation Survey",
                            "keywords": "Kuiper Belt, Trans Neptunian Object, Occultations"
                            },
              "target": {"name": f"{fobj['header']['pointing']['field_id']}",
                         "type": f"{fobj['header']['run']['imgtype']}",
                         "targetID": f"{fobj['haeder']['object']['obj_type']}:{fobj['header']['object']['obj_id']}",
                         "type": f"{fobj['header']['object']['obj_type']}"
                         },
              # "telescope": {"name": f"{fobj['header']['device']}"
              "telescope": {
                  "name": "TAOS-II",
                  # Long = -115.454
                  # Lat  =   31.041
                  # Elev = 2820m
                  "geoLocationX": "",
                  "geoLocationY": "",
                  "geoLocationZ": "",
                        }
              "environment": {
                  "ambientTemp": "db.taos2.temperature"
              }
          "Plane": {
              "productID": ""
          }
        }

        :param bp:
        :return:
        """
        """
        need n WCS instances - one per site, so three per file, so, how
        to know when to create those, and how many to create?
        - I think that knowing there are n (three) per file is just a thing 
        that is known ahead of time, and that the path to the 'n' bits is
        part of the construction of HCF5Parser
        - then, the separate WCSParser bits are constructed as part of that
        handling?
        
        For FITS:
        e.g. Chunk.position.axis.axis1.ctype = CTYPE1

        For HDF5:
        e.g. Chunk.position.axis.axis1.ctype = /header/wcs/ctype[0]
        
        BUT
        - need an astropy.wcs.WCS construction for correctness/consistency

        - for the FITS file, the WCS construction is handled by astropy
        - for the HDF5 file, the WCS construction needs to be handled by
          blueprint (?) code, so something like:
          w.wcs.ctype = [value from HDF5 file, which is found by looking 
                         up the CAOM2 key, which is then used to retrieve
                         the value from the HDF5 file]
        """

        super().accumulate_blueprint(bp)
        bp.configure_time_axis(3)
        bp.configure_energy_axis(4)

        self._prefix = '//header'
        if self._storage_name.is_lightcurve:
            self._prefix = '//obs/header'
        utc_now = datetime.utcnow()
        bp.set('Observation.intent', self.get_observation_intent(0))
        bp.set('Plane.dataProductType', DataProductType.IMAGE)

        bp.set_default('Observation.metaRelease', utc_now)
        bp.set('Plane.dataRelease', utc_now)
        bp.set('Plane.metaRelease', utc_now)

        bp.set('Observation.type', ([f'{self._prefix}/run/obstype'], 'IMAGE'))
        bp.set('Observation.instrument.name', 'get_instrument_name()')
        bp.set('Observation.instrument.keywords', ([f'{self._prefix}/run/exptype'], None))
        bp.set('Observation.proposal.id', ([f'{self._prefix}/run/run_seq'], None))
        bp.set('Observation.proposal.pi', ([f'{self._prefix}/run/observer'], None))
        bp.set('Observation.proposal.project', 'TAOS2')
        bp.set('Observation.proposal.title', 'Transneptunian Automated Occultation Survey')
        bp.set('Observation.proposal.keywords', set(['Kuiper Belt', 'Trans Neptunian Object', 'Occultations']))

        bp.set('Observation.telescope.name', self.get_telescope_name(0))
        # x, y, z = ac.get_location(31.041, -115.454, 2820)
        bp.set('Observation.telescope.geoLocationX', -2351818.5502637075)
        bp.set('Observation.telescope.geoLocationY', -4940894.8697881885)
        bp.set('Observation.telescope.geoLocationZ', 3271243.2086214763)

        if 'master' in self._storage_name.file_name:
            bp.set('Observation.algorithm.name', 'master')
            bp.set('DerivedObservation.members', {})
            bp.set('Plane.provenance.inputs', ([f'{self._prefix}/input_files'], None))
            bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)

        bp.set('Artifact.productType', self._storage_name.get_product_type)

        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.exposure', ([f'{self._prefix}/exposure/exposure'], None))
        bp.set('Chunk.time.mjdref', ([f'{self._prefix}/coord_params/mjdref'], None))
        # JJK 30-01-23
        # timesys is UTC.
        # bp.set('Chunk.time.timesys', 'UTC')

        # SGw - 23-07-22
        # Detector energy information from Figure 4 here:
        # https://www.spiedigitallibrary.org/proceedings/Download?urlId=10.1117%2F12.2561204
        #
        # This reference says there is no filter:
        # https://www.spiedigitallibrary.org/conference-proceedings-of-spie/9908/1/
        # The-prototype-cameras-for-trans-Neptunian-automatic-occultation-survey/10.1117/12.2232062.full
        #
        # FWHM => 430nm to 830nm
        #
        # JJK 30-01-23
        # resolving power: R == Lambda/Delta_Lambda
        # Where Lambda is the wavelength at the middle of the bandpass and Delta_Lambda is width.  In this case, Lambda
        # = (8.3E-7 + 4.3E-7)/2 and Delta_lambda = 8.3E-7 - 4.3E-7
        # So+   (8.3+4.3)/(2*(8.3-4.3)) = 1.575
        #
        bp.set('Chunk.energy.specsys', 'TOPOCENT')
        bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
        bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
        bp.set('Chunk.energy.bandpassName', 'CLEAR')
        bp.set('Chunk.energy.resolvingPower', 1.575)
        bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
        # units as output by astropy.wcs.WCS
        bp.set('Chunk.energy.axis.axis.cunit', 'm')

    def get_instrument_name(self, ext):
        if self._storage_name.is_fsc:
            result = 'TAOSZWOCAM'
        else:
            result = self._lookup(f'{self._prefix}/run/origin')
        return result

    def get_observation_intent(self, ext):
        result = ObservationIntentType.CALIBRATION
        if self._storage_name.get_product_type == ProductType.SCIENCE:
            result = ObservationIntentType.SCIENCE
        return result

    def get_telescope_name(self, ext):
        result = 'TAOSII'
        if self._storage_name.is_fsc:
            result = 'TAOSIIFSC'
        return result

    def _lookup(self, things):
        things_replaced = things.replace('//', '', 1)
        def hd5f_visit(name, object):
            if (
                isinstance(object, h5py.Dataset)
                and object.dtype.names is not None
            ):
                # 'name' starts without a directory separator
                if things_replaced.startswith(name):
                    for d_name in object.dtype.names:
                        temp = f'{name}/{d_name}'
                        if temp == things_replaced:
                            return object[d_name]

        result = self._h5_file.visititems(hd5f_visit)
        if result is None:
            self._logger.warning(f'Could not find {things} in {self._storage_name.file_name}')
        else:
            self._logger.debug(f'Found {result} for {things}')
            if hasattr(result, 'decode'):
                result = result.decode('utf-8')
        return result

    def _update_artifact(self, artifact):
        self._logger.debug('Begin _update_artifact')
        for part in artifact.parts.values():
            for chunk in part.chunks:
                if (
                    chunk.time is not None
                    and chunk.time.axis is not None
                    and chunk.time.axis.function is not None
                ):
                    if chunk.time.axis.range is not None:
                        # time is range, not function, and at this time, the
                        # blueprint will always set the function
                        chunk.time.axis.function = None
                    # chunk.time.trefpos = None
                    chunk.time_axis = None
                if (
                    chunk.energy is not None
                    and chunk.energy.axis is not None
                    and chunk.energy.axis.function is not None
                ):
                    if chunk.energy.axis.range is not None:
                        # energy is range, not function, and at this time, the
                        # blueprint will always set the function
                        chunk.energy.axis.function = None
                    chunk.energy_axis = None
                if (
                    chunk.position is not None
                    and chunk.position_axis_1 is not None
                ):
                    chunk.position_axis_1 = None
                    chunk.position_axis_2 = None
                if self._storage_name.is_lightcurve and chunk.position is not None:
                    chunk.position = None
        self._logger.debug('End _update_artifact')


class DomeflatMapping(BasicMapping):

    def __init__(self, storage_name, h5file, clients, observable, observation):
        super().__init__(storage_name, h5file, clients, observable, observation)

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Chunk.time.axis.range.start.pix', 0.5)
        bp.set('Chunk.time.axis.range.end.pix', 1.5)
        bp.set('Chunk.time.axis.range.start.val', ([f'{self._prefix}/exposure/mjdstart'], None))
        bp.set('Chunk.time.axis.range.end.val', ([f'{self._prefix}/exposure/mjdend'], None))

        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        bp.set('Chunk.energy.axis.range.start.val', 4.3e-7)
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        bp.set('Chunk.energy.axis.range.end.val', 8.3e-7)


class SingleMapping(BasicMapping):

    def __init__(self, storage_name, h5file, clients, observable, observation):
        super().__init__(storage_name, h5file, clients=clients, observable=observable, observation=observation)

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)

        bp.set('Observation.telescope.geoLocationX', ([f'{self._prefix}/avg_location/obsgeo(0)'], None))
        bp.set('Observation.telescope.geoLocationY', ([f'{self._prefix}/avg_location/obsgeo(1)'], None))
        bp.set('Observation.telescope.geoLocationZ', ([f'{self._prefix}/avg_location/obsgeo(2)'], None))
        bp.set('Observation.target.name', ([f'{self._prefix}/pointing/field_id'], None))

        bp.set('Artifact.productType', self._storage_name.get_product_type)

        bp.set('Chunk.time.axis.function.refCoord.pix', (['/header/wcs/crpix(2)'], None))
        bp.set('Chunk.time.axis.function.refCoord.val', (['/header/wcs/crval(2)'], None))
        bp.set('Chunk.time.axis.axis.ctype', (['/header/wcs/ctype(2)'], None))
        bp.set('Chunk.time.axis.axis.cunit', (['/header/wcs/cunit(2)'], None))
        bp.set('Chunk.time.timesys', (['/header/wcs/cname(2)'], None))

        bp.set('Chunk.energy.axis.function.refCoord.pix', (['/header/wcs/crpix(3)'], None))
        bp.set('Chunk.energy.axis.function.refCoord.val', (['/header/wcs/crval(3)'], None))
        bp.set('Chunk.energy.axis.axis.ctype', (['/header/wcs/ctype(3)'], None))
        bp.set('Chunk.energy.axis.axis.cunit', (['/header/wcs/cunit(3)'], None))


class PointingMapping(SingleMapping):

    def __init__(self, storage_name, h5file, clients, observable, observation):
        super().__init__(storage_name, h5file, clients, observable, observation)

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((1, 2))
        bp.set('Plane.calibrationLevel', CalibrationLevel.RAW_STANDARD)

        bp.set('Chunk.position.axis.function.dimension.naxis1', 1920)
        bp.set('Chunk.position.axis.function.dimension.naxis2', 4608)
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix', (['/header/wcs/crpix(0)'], None))
        bp.set('Chunk.position.axis.function.refCoord.coord1.val', (['/header/wcs/crval(0)'], None))
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix', (['/header/wcs/crpix(1)'], None))
        bp.set('Chunk.position.axis.function.refCoord.coord2.val', (['/header/wcs/crval(1)'], None))
        bp.set('Chunk.position.axis.axis1.ctype', (['/header/wcs/ctype(0)'], None))
        bp.set('Chunk.position.axis.axis1.cunit', (['/header/wcs/cunit(0)'], None))
        bp.set('Chunk.position.axis.axis2.ctype', (['/header/wcs/ctype(1)'], None))
        bp.set('Chunk.position.axis.axis2.cunit', (['/header/wcs/cunit(1)'], None))
        bp.set('Chunk.position.axis.function.cd11', (['/header/wcs/cd(0:0)'], None))
        bp.set('Chunk.position.axis.function.cd12', (['/header/wcs/cd(0:1)'], None))
        bp.set('Chunk.position.axis.function.cd21', (['/header/wcs/cd(1:0)'], None))
        bp.set('Chunk.position.axis.function.cd22', (['/header/wcs/cd(1:1)'], None))
        bp.set('Chunk.position.equinox', ([f'{self._prefix}/object/equinox'], None))
        # JJK 30-01-23
        # coordsys: null should likely be coordsys: ICRS
        bp.set('Chunk.position.coordsys', ([f'{self._prefix}/coord_params/radecsys'], None))
        # logging.error(bp)


class TimeseriesMapping(PointingMapping):

    def __init__(self, storage_name, h5file, clients, observable, observation):
        super().__init__(storage_name, h5file, clients, observable, observation)

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Plane.calibrationLevel', CalibrationLevel.PRODUCT)
        bp.set('Plane.dataProductType', DataProductType.TIMESERIES)

        bp.set('Observation.target.targetID', ([f'{self._prefix}/object/obj_id'], None))
        bp.set('Observation.target.type', ([f'{self._prefix}/obj/obj_type'], None))
        bp.set('Observation.target_position.point.cval1', 'get_target_position_cval1()')
        bp.set('Observation.target_position.point.cval2', 'get_target_position_cval2()')
        bp.set('Observation.target_position.coordsys', 'FK5')
        bp.set('Observation.target_position.equinox', ([f'{self._prefix}/object/equinox'], None))

    def _get_target_position(self):
        obj_ra = self._lookup(f'{self._prefix}/object/obj_ra')
        obj_dec = self._lookup(f'{self._prefix}/object/obj_dec')
        if obj_ra is None or obj_dec is None:
            ra = None
            dec = None
        else:
            ra, dec = build_ra_dec_as_deg(obj_ra, obj_dec)
        return ra, dec

    def get_target_position_cval1(self, ext):
        ra, dec_ignore = self._get_target_position()
        return ra

    def get_target_position_cval2(self, ext):
        ra_ignore, dec = self._get_target_position()
        return dec


class TAOSIIName(StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage
    """

    TAOSII_NAME_PATTERN = '*'

    def __init__(self, file_name=None, source_names=None):

        super().__init__(file_name=basename(file_name), source_names=source_names)

    @property
    def get_product_type(self):
        if '_star' in self._file_name or '_fsc' in self._file_name:
            return ProductType.SCIENCE
        else:
            return ProductType.CALIBRATION

    @property
    def is_fsc(self):
        return '_fsc' in self._file_name

    @property
    def is_lightcurve(self):
        return self._product_id.startswith('taos2lcv_')

    def is_valid(self):
        return True

    def set_file_id(self):
        self._file_id = StorageName.remove_extensions(self._file_name).replace('.hdf5', '').replace('.h5', '')

    def set_obs_id(self):
        self._obs_id = self._file_id.replace('lcv_obs', '')


class TAOSII2caom2Visitor(Fits2caom2Visitor):

    def __init__(self, observation, **kwargs):
        super().__init__(observation, **kwargs)
        self._h5_file = h5py.File(self._storage_name.source_names[0])

    def _get_blueprint(self, instantiated_class):
        return Hdf5ObsBlueprint(instantiated_class=instantiated_class)

    def _get_mapping(self, headers):
        if '_star' in self._storage_name.file_name:
            result = TimeseriesMapping(
                self._storage_name, self._h5_file, self._clients, self._observable, self._observation
            )
        elif (
            '_focus' in self._storage_name.file_name
            or '_point' in self._storage_name.file_name
            or self._storage_name.is_fsc
        ):
            result = PointingMapping(
                self._storage_name, self._h5_file, self._clients, self._observable, self._observation
            )
        elif (
            '_domeflat' in self._storage_name.file_name
            or '_dark' in self._storage_name.file_name
            or '_bias' in self._storage_name.file_name
        ):
            result = DomeflatMapping(
                self._storage_name, self._h5_file, self._clients, self._observable, self._observation
            )
        else:
            result = SingleMapping(
                self._storage_name, self._h5_file, self._clients, self._observable, self._observation
            )
        self._logger.debug(f'Created mapping {result.__class__.__name__}.')
        return result

    def _get_parser(self, headers, blueprint, uri):
        self._logger.debug(
            f'Using an Hdf5Parser for {self._storage_name.file_uri} '
        )
        # this assumes only working with one file at a time, and also, that
        # the file is local (which, as of the time of this writing, the file
        # has to be local, until there is an --fhead option for HDF5 files)
        return Hdf5Parser(blueprint, uri, self._h5_file, 'images')

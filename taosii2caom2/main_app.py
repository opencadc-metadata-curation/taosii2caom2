# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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


Mapping Notes:
1. Single Lightcurve File - up to three lightcurves, one from each telescope.
    a. EXPTYPE - generally used for algorithm.name, not expected to be found.
               - algorithm.name is hard-coded
    b. DATASEC - generally used for finding NAXISij values, not expected to be found.
               - refCoord.coord1.pix, refCoord.coord2.pix hard-coded to 0.5
2. High-speed Window Mode

ML - Matt Lehner
JJK - JJ Kavelaars

"""

import h5py
import logging

from collections import defaultdict
from datetime import datetime
from dateutil.relativedelta import relativedelta
from os.path import basename, join

from astropy.time import Time

from caom2 import CalibrationLevel, DataProductType, ProductType
from caom2 import ObservationIntentType, ObservationURI, PlaneURI, TypedSet

from caom2utils.caom2blueprint import Hdf5ObsBlueprint, Hdf5Parser
from caom2utils.blueprints import _to_float
from caom2utils.data_util import get_local_file_info
from caom2pipe.astro_composable import build_ra_dec_as_deg
from caom2pipe.caom_composable import Fits2caom2VisitorRunnerMeta, TelescopeMapping2
from caom2pipe.execute_composable import (
    CaomExecuteRunnerMeta,
    OrganizeExecutesRunnerMeta,
    NoFheadScrapeRunnerMeta,
    NoFheadStoreVisitRunnerMeta,
    NoFheadVisitRunnerMeta,
)
from caom2pipe.manage_composable import CadcException, CaomName, make_datetime, StorageName, TaskType, ValueRepairCache


__all__ = ['set_storage_name_from_local_preconditions', 'TAOSII2caom2Visitor', 'TAOSIIName']


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


class NoWcsMapping(TelescopeMapping2):

    def __init__(self, storage_name, clients, reporter, observation, config, prefix):
        super().__init__(
            storage_name, clients=clients, reporter=reporter, observation=observation, config=config
        )
        self._prefix = f'//{prefix}'

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
        bp.set('Observation.intent', self.get_observation_intent(0))
        bp.set_default('Observation.metaRelease', self._get_meta_release_value())
        bp.set('Observation.sequenceNumber', ([f'{self._prefix}/HEADER/RUN/RUN_SEQ'], None))
        bp.set('Observation.type', ([f'{self._prefix}/OBS/OBSTYPE'], 'IMAGE'))
        bp.set('Observation.instrument.name', 'get_instrument_name()')
        bp.set('Observation.instrument.keywords', ([f'{self._prefix}/HEADER/OBS/EXPTYPE'], None))
        bp.set('Observation.proposal.id', ([f'{self._prefix}/HEADER/RUN/RUN_SEQ'], None))
        bp.set('Observation.proposal.pi', ([f'{self._prefix}/HEADER/RUN/OBSERVER'], None))
        bp.set('Observation.proposal.project', 'TAOS2')
        bp.set('Observation.proposal.title', 'Transneptunian Automated Occultation Survey')
        bp.set('Observation.proposal.keywords', set(['Kuiper Belt', 'Trans Neptunian Object', 'Occultations']))
        bp.set('Observation.telescope.name', self.get_telescope_name(0))
        # x, y, z = ac.get_location(31.041, -115.454, 2820)
        bp.set('Observation.telescope.geoLocationX', -2351818.5502637075)
        bp.set('Observation.telescope.geoLocationY', -4940894.8697881885)
        bp.set('Observation.telescope.geoLocationZ', 3271243.2086214763)
        # all the samples are calibrated
        bp.set('Observation.algorithm.name', self.get_algorithm_name())
        bp.set('DerivedObservation.members', {})

        bp.set('Plane.dataRelease', self._get_data_release())
        bp.set('Plane.metaRelease', self._get_meta_release_value())
        bp.set('Plane.dataProductType', DataProductType.TIMESERIES)
        bp.set('Plane.provenance.name', (['//FILE/ORIGIN'], None))
        bp.set('Plane.provenance.version', self.get_provenance_version())
        bp.set('Plane.provenance.lastExecuted', (['//FILE/FILE_DATE'], None))
        bp.set('Plane.provenance.inputs', self.get_provenance_inputs())
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)

        bp.set('Artifact.productType', self.get_product_type())

    def _get_data_release(self):
        result = None
        meta_release_value = self._get_meta_release_value()
        if meta_release_value:
            # JJK, ML 23-10-24 - 5 year proprietary period for the data
            # JJK 29-10-24 - no metadata proprietary period
            x = relativedelta(years=5)
            y = datetime.fromisoformat(meta_release_value)
            result = x + y
        return result

    def _get_meta_release(self):
        result = None
        mjdref = self._lookup_for_links('HEADER/COORDSYS/MJDREF')
        if mjdref:
            obs_start = _to_float(self._lookup_for_links('HEADER/EXPOSURE/OBSSTART'))
            if obs_start:
                result = Time(mjdref, obs_start, format='mjd', scale='tt')
        return result

    def _get_meta_release_value(self):
        result = self._get_meta_release()
        if result:
            result.format = 'isot'
            result = result.value
        return result

    def _get_obs_type(self):
        return self._lookup(f'{self._prefix}/HEADER/OBS/OBSTYPE')

    def _get_provenance_inputs(self):
        return self._lookup(f'{self._prefix}/HEADER/CALIBRATION/FILE')

    def get_algorithm_name(self):
        result = 'lightcurve'
        temp = self._lookup(f'{self._prefix}/HEADER/OBS/EXPTYPE')
        if temp:
            result = temp.lower()
        return result

    def get_instrument_name(self, ext):
        if self._storage_name.is_fsc:
            result = 'TAOSZWOCAM'
        else:
            result = self._lookup(f'//FILE/ORIGIN')
        return result

    def get_observation_intent(self, ext):
        if self._storage_name.is_lightcurve:
            result = ObservationIntentType.SCIENCE
        else:
            obstype = self._get_obs_type()
            result = ObservationIntentType.CALIBRATION
            if obstype == 'IMAGE':
                result = ObservationIntentType.SCIENCE
        return result

    def get_product_type(self):
        if self._storage_name.is_lightcurve:
            result = ProductType.SCIENCE
        else:
            obstype = self._get_obs_type()
            x = {
                'IMAGE': ProductType.SCIENCE,
                'FOCUS': ProductType.CALIBRATION,
                'POINTING': ProductType.CALIBRATION,
                'POINTING_MODEL': ProductType.CALIBRATION,
                'BIAS': ProductType.BIAS,
                'DARK': ProductType.DARK,
                # avoid server error:
                # "ERROR: invalid input: ... reason: XML failed schema validation: Error on line 43:
                # cvc-enumeration-valid: Value 'flat' is not facet-valid with respect to enumeration
                # '[science, calibration, auxiliary, info, preview, catalog, noise, weight, thumbnail]'.
                # It must be a value from the enumeration."
                'DOMEFLAT': ProductType.CALIBRATION,
                'SKYFLAT': ProductType.CALIBRATION,
            }
            result = x.get(obstype, ProductType.AUXILIARY)
        return result

    def get_provenance_inputs(self):
        temp = self._get_provenance_inputs()
        plane_inputs = TypedSet(PlaneURI)
        if temp is not None:
            for entry in temp:
                obs_id = TAOSIIName.replace_for_obs_id(entry.decode('utf-8'))
                obs_member_uri_str = CaomName.make_obs_uri_from_obs_id(self._storage_name.collection, obs_id)
                obs_member_uri = ObservationURI(obs_member_uri_str)
                plane_uri = PlaneURI.get_plane_uri(obs_member_uri, entry.decode('utf-8'))
                plane_inputs.add(plane_uri)
                self._logger.debug(f'Adding PlaneURI {plane_uri}')
        return plane_inputs

    def get_provenance_version(self):
        temp = self._lookup('//FILE/VERSION')
        return '.'.join(str(ii) for ii in temp)

    def get_telescope_name(self, ext):
        result = 'TAOSII'
        if self._storage_name.is_fsc:
            result = 'TAOSIIFSC'
        return result

    def _lookup(self, things, fishing=False):
        """
        :params
        :things str - POSIX path breadcrumbs
        :fishing bool - True if doing a lookup for characterization, so set the logging level accordingly
        """
        result = _lookup(self._storage_name.descriptors.get(self._storage_name.file_uri), things)
        if result is None:
            msg = f'Could not find {things} in {self._storage_name.file_name}'
            if fishing:
                self._logger.debug(msg)
            else:
                self._logger.warning(msg)
        else:
            self._logger.debug(f'Found {result} for {things}')
        return result

    def _lookup_for_links(self, key_suffix):
        result = self._lookup(f'{self._prefix}/{key_suffix}')
        if result is None:
            result = self._lookup(f'LIGHTCURVE/{key_suffix}')
        return result

    def _update_artifact(self, artifact):
        self._logger.debug('Begin _update_artifact')
        for part in artifact.parts.values():
            for chunk in part.chunks:
                if chunk.time is not None and chunk.time.axis is not None and chunk.time.axis.function is not None:
                    chunk.time_axis = None
                if (
                    chunk.energy is not None
                    and chunk.energy.axis is not None
                    and chunk.energy.axis.function is not None
                ):
                    chunk.energy_axis = None
                if chunk.position is not None and chunk.position_axis_1 is not None:
                    chunk.position_axis_1 = None
                    chunk.position_axis_2 = None
        self._logger.debug('End _update_artifact')


class BasicMapping(NoWcsMapping):

    def __init__(
        self, storage_name, clients, reporter, observation, config, prefix, time_axis=3, energy_axis=4
    ):
        super().__init__(storage_name, clients, reporter, observation, config, prefix)
        self._time_axis = time_axis
        self._energy_axis = energy_axis

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.configure_time_axis(self._time_axis)
        bp.configure_energy_axis(self._energy_axis)

        bp.set('Chunk.time.axis.axis.ctype', ([f'/HEADER/WCS/CTYPE({self._time_axis - 1})'], None))
        bp.set('Chunk.time.axis.axis.cunit', ([f'/HEADER/WCS/CUNIT({self._time_axis - 1})'], None))
        bp.set(
            'Chunk.time.axis.function.delta', ([f'/HEADER/WCS/CD({self._time_axis - 1}:{self._time_axis - 1})'], None)
        )
        bp.set('Chunk.time.axis.function.refCoord.pix', ([f'/HEADER/WCS/CRPIX({self._time_axis - 1})'], None))
        bp.set('Chunk.time.axis.function.refCoord.val', ([f'/HEADER/WCS/CRVAL({self._time_axis - 1})'], None))
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', ([f'/HEADER/WCS/CUNIT({self._time_axis - 1})'], None))
        bp.set('Chunk.time.timesys', ([f'/HEADER/WCS/CTYPE({self._time_axis - 1})'], None))
        bp.set('Chunk.time.exposure', ([f'{self._prefix}/HEADER/EXPOSURE/EXPOSURE'], None))
        bp.set('Chunk.time.mjdref', ([f'{self._prefix}/HEADER/COORDSYS/MJDREF'], None))
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
        bp.set('Chunk.energy.axis.axis.ctype', ([f'/HEADER/WCS/CTYPE({self._energy_axis - 1})'], None))
        bp.set('Chunk.energy.axis.axis.cunit', ([f'/HEADER/WCS/CUNIT({self._energy_axis - 1})'], None))
        bp.set(
            'Chunk.energy.axis.function.delta',
            ([f'/HEADER/WCS/CD({self._energy_axis - 1}:{self._energy_axis - 1})'], None),
        )
        bp.set('Chunk.energy.axis.function.refCoord.pix', ([f'/HEADER/WCS/CRPIX({self._energy_axis - 1})'], None))
        bp.set('Chunk.energy.axis.function.refCoord.val', ([f'/HEADER/WCS/CRVAL({self._energy_axis - 1})'], None))
        bp.set('Chunk.energy.axis.axis.ctype', ([f'/HEADER/WCS/CTYPE({self._energy_axis - 1})'], None))
        bp.set('Chunk.energy.axis.axis.cunit', ([f'/HEADER/WCS/CUNIT({self._energy_axis - 1})'], None))
        bp.set('Chunk.energy.specsys', 'TOPOCENT')
        bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
        bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
        bp.set('Chunk.energy.bandpassName', 'CLEAR')
        bp.set('Chunk.energy.resolvingPower', 1.575)

    def _get_time_exposure(self):
        return _to_float(self._lookup_for_links('HEADER/EXPOSURE/EXPOSURE'))


class DomeflatMapping(NoWcsMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)

        bp.configure_time_axis(1)
        bp.configure_energy_axis(2)

        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.range.start.pix', 0.5)
        bp.set('Chunk.time.axis.range.end.pix', 1.5)
        bp.set('Chunk.time.axis.range.start.val', ([f'{self._prefix}/HEADER/EXPOSURE/OBSSTART'], None))
        bp.set('Chunk.time.axis.range.end.val', ([f'{self._prefix}/HEADER/EXPOSURE/OBSEND'], None))
        bp.set('Chunk.time.timesys', ([f'{self._prefix}/HEADER/COORDSYS/TIMESYS'], None))
        bp.set('Chunk.time.exposure', ([f'{self._prefix}/HEADER/EXPOSURE/EXPOSURE'], None))
        bp.set('Chunk.time.mjdref', ([f'{self._prefix}/HEADER/COORDSYS/MJDREF'], None))

        bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
        bp.set('Chunk.energy.axis.axis.cunit', 'm')
        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        # values from ML - 23-10-24
        bp.set('Chunk.energy.axis.range.start.val', 4.0e-7)
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        bp.set('Chunk.energy.axis.range.end.val', 7.2e-7)
        bp.set('Chunk.energy.specsys', 'TOPOCENT')
        bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
        bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
        bp.set('Chunk.energy.bandpassName', 'CLEAR')
        bp.set('Chunk.energy.resolvingPower', 1.575)


class Pointing(BasicMapping):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Observation.telescope.geoLocationX', ([f'{self._prefix}/HEADER/LOCATION/OBSGEO(0)'], None))
        bp.set('Observation.telescope.geoLocationY', ([f'{self._prefix}/HEADER/LOCATION/OBSGEO(1)'], None))
        bp.set('Observation.telescope.geoLocationZ', ([f'{self._prefix}/HEADER/LOCATION/OBSGEO(2)'], None))
        bp.set('Observation.target.name', ([f'{self._prefix}/HEADER/POINTING/FIELD_ID'], None))

        bp.set('Artifact.productType', self.get_product_type())
        self._logger.debug('End accumulate_blueprint')


class TargetSpectralTemporal(Pointing):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Observation.target.targetID', self.get_target_id())
        bp.set('Observation.target.type', self.get_target_type())
        bp.set('Observation.target_position.point.cval1', self.get_target_position_cval1())
        bp.set('Observation.target_position.point.cval2', self.get_target_position_cval2())
        bp.set('Observation.target_position.coordsys', self.get_target_coordsys())

    def _get_target_position(self):
        obj_ra = self._lookup(f'{self._prefix}/HEADER/OBJECT/OBJ_RA')
        obj_dec = self._lookup(f'{self._prefix}/HEADER/OBJECT/OBJ_DEC')
        if obj_ra is None or obj_dec is None:
            ra = None
            dec = None
        else:
            ra, dec = build_ra_dec_as_deg(obj_ra, obj_dec)
        return ra, dec

    def get_target_coordsys(self):
        return self._lookup(f'{self._prefix}/HEADER/COORDSYS/RADECSYS')

    def get_target_id(self):
        return self._lookup(f'{self._prefix}/HEADER/OBJECT/OBJ_ID')

    def get_target_position_cval1(self):
        ra, _ = self._get_target_position()
        return ra

    def get_target_position_cval2(self):
        _, dec = self._get_target_position()
        return dec

    def get_target_type(self):
        obj_type = self._lookup(f'{self._prefix}/HEADER/OBJECT/OBJTYPE')
        result = 'object'
        if obj_type != 'TAOS_FIELDSTAR':
            self._logger.warning(f'Unexpected OBJTYPE {obj_type} for {self._storage_name.file_uri}')
            result = None
        return result


class Single(TargetSpectralTemporal):

    def __init__(
        self,
        storage_name,
        clients,
        reporter,
        observation,
        config,
        prefix,
        time_axis,
        energy_axis,
        spatial_axis_1,
    ):
        super().__init__(
            storage_name, clients, reporter, observation, config, prefix, time_axis, energy_axis
        )
        # expect this to be 1 for 0, 1 CD matrix references
        #                   2 for 1, 2 CD matrix references
        self._spatial_axis = spatial_axis_1

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.configure_position_axes((self._spatial_axis - 1, self._spatial_axis))
        bp.set('Plane.calibrationLevel', CalibrationLevel.RAW_STANDARD)

        bp.set('Chunk.position.axis.function.dimension.naxis1', self.get_naxis1())
        bp.set('Chunk.position.axis.function.dimension.naxis2', self.get_naxis2())
        bp.set(
            'Chunk.position.axis.function.refCoord.coord1.pix',
            ([f'/HEADER/WCS/CRPIX({self._spatial_axis - 1})'], None),
        )
        bp.set(
            'Chunk.position.axis.function.refCoord.coord1.val',
            ([f'/HEADER/WCS/CRVAL({self._spatial_axis - 1})'], None),
        )
        bp.set(
            'Chunk.position.axis.function.refCoord.coord2.pix',
            ([f'/HEADER/WCS/CRPIX({self._spatial_axis})'], None),
        )
        bp.set(
            'Chunk.position.axis.function.refCoord.coord2.val',
            ([f'/HEADER/WCS/CRVAL({self._spatial_axis})'], None),
        )
        bp.set('Chunk.position.axis.axis1.ctype', ([f'/HEADER/WCS/CTYPE({self._spatial_axis - 1})'], None))
        bp.set('Chunk.position.axis.axis1.cunit', ([f'/HEADER/WCS/CUNIT({self._spatial_axis - 1})'], None))
        bp.set('Chunk.position.axis.axis2.ctype', ([f'/HEADER/WCS/CTYPE({self._spatial_axis})'], None))
        bp.set('Chunk.position.axis.axis2.cunit', ([f'/HEADER/WCS/CUNIT({self._spatial_axis})'], None))
        bp.set(
            'Chunk.position.axis.function.cd11',
            ([f'/HEADER/WCS/CD({self._spatial_axis - 1}:{self._spatial_axis - 1})'], None),
        )
        bp.set(
            'Chunk.position.axis.function.cd12',
            ([f'/HEADER/WCS/CD({self._spatial_axis - 1}:{self._spatial_axis})'], None),
        )
        bp.set(
            'Chunk.position.axis.function.cd21',
            ([f'/HEADER/WCS/CD({self._spatial_axis}:{self._spatial_axis - 1})'], None),
        )
        bp.set(
            'Chunk.position.axis.function.cd22',
            ([f'/HEADER/WCS/CD({self._spatial_axis}:{self._spatial_axis})'], None),
        )
        # JJK 30-01-23
        # coordsys: null should likely be coordsys: ICRS
        bp.set('Chunk.position.coordsys', ([f'{self._prefix}/HEADER/COORDSYS/RADECSYS'], None))

    def _get_naxis(self, index):
        result = None
        temp = self._lookup(f'{self._prefix}/SITE1/HEADER/PIXSEC/DATASEC')
        if not temp:
            temp = self._lookup(f'{self._prefix}/SITE_1/HEADER/PIXSEC/DATASEC')
        if temp:
            result = temp[index][1] - temp[index][0] + 1
        return result

    def get_naxis1(self):
        return self._get_naxis(0)

    def get_naxis2(self):
        return self._get_naxis(1)


class FullFrameImage(Single):
    def __init__(self, storage_name, clients, reporter, observation, config, prefix):
        super().__init__(
            storage_name,
            clients,
            reporter,
            observation,
            config,
            prefix,
            time_axis=3,
            energy_axis=4,
            spatial_axis_1=1,
        )

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        bp.set('Plane.dataProductType', DataProductType.IMAGE)

    def _get_target_position(self):
        pointing_ra = self._lookup(f'{self._prefix}/HEADER/POINTING/TEL_RA')
        pointing_dec = self._lookup(f'{self._prefix}/HEADER/POINTING/TEL_DEC')
        if pointing_ra is None or pointing_dec is None:
            ra = None
            dec = None
        else:
            ra, dec = build_ra_dec_as_deg(pointing_ra, pointing_dec)
        return ra, dec

    def get_target_coordsys(self):
        return self._lookup(f'{self._prefix}/HEADER/COORDSYS/RADECSYS')

    def get_target_id(self):
        return self._lookup(f'{self._prefix}/HEADER/POINTING/FIELD')

    def get_target_position_cval1(self):
        ra, _ = self._get_target_position()
        return ra

    def get_target_position_cval2(self):
        _, dec = self._get_target_position()
        return dec

    def get_target_type(self):
        return 'field'


class TimeseriesMapping(Single):

    def __init__(self, storage_name, clients, reporter, observation, config, prefix):
        super().__init__(
            storage_name,
            clients,
            reporter,
            observation,
            config,
            prefix,
            time_axis=3,
            energy_axis=4,
            spatial_axis_1=1,
        )

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.set('Plane.calibrationLevel', CalibrationLevel.PRODUCT)


class Lightcurve(Single):
    def __init__(
        self,
        storage_name,
        clients,
        reporter,
        observation,
        config,
        prefix,
        time_axis=1,
        energy_axis=4,
        spatial_axis_1=2,
    ):
        super().__init__(
            storage_name,
            clients,
            reporter,
            observation,
            config,
            prefix,
            time_axis,
            energy_axis,
            spatial_axis_1,
        )

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        # guess based on other test light-curve file cal levels
        bp.set('Plane.calibrationLevel', CalibrationLevel.ANALYSIS_PRODUCT)
        bp.set('Plane.dataProductType', DataProductType.TIMESERIES)
        bp.set('Artifact.productType', ProductType.SCIENCE)

        # this is my guess
        bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
        bp.set('Chunk.position.axis.function.dimension.naxis2', 1)
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix', 0.0)
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 0.0)

        bp.set('Chunk.time.axis.function.naxis', self._get_lightcurve_time_naxis())
        bp.set('Chunk.time.exposure', self._get_time_exposure())
        bp.set('Chunk.time.mjdref', self._get_time_mjdref())

    def _get_lightcurve_time_naxis(self):
        # set the value as a function to avoid relying on the Part index
        result = self._lookup(f'{self._prefix}/HEADER/EXPOSURE/NUM_EPOCHS')
        if result is None:
            result = self._lookup(f'IMAGE/HEADER/EXPOSURE/NUM_EPOCHS')
        return result

    def _get_provenance_inputs(self):
        return self._lookup(f'{self._prefix}/HEADER/CAL/FILE')

    def _get_time_exposure(self):
        # the default option is for the stand-alone lightcurve files
        result = self._lookup(f'{self._prefix}/HEADER/EXPOSURE/EXPOSURE')
        if result is None:
            # this option is for the files that have images and lightcurves together,
            # and the HEADER metadata is a link
            result = self._lookup(f'IMAGE/HEADER/EXPOSURE/EXPOSURE')
        return result

    def _get_time_mjdref(self):
        # the default option is for the stand-alone lightcurve files
        result = self._lookup(f'{self._prefix}/HEADER/COORDSYS/MJDREF')
        if result is None:
            # this option is for the files that have images and lightcurves together,
            # and the HEADER metadata is a link
            result = self._lookup(f'IMAGE/HEADER/COORDSYS/MJDREF')
        return result


class FSC(FullFrameImage):

    def accumulate_blueprint(self, bp):
        super().accumulate_blueprint(bp)
        bp.clear('Chunk.energy.resolvingPower')


class TAOSIIName(StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage
    """

    TAOSII_NAME_PATTERN = '*'

    def __init__(self, source_names=None):
        # the descriptors values get assigned in the _set_preconditions call the execute_composable.CaomExecutes
        # specialization
        self._file_uri = None
        self.descriptors = {}
        super().__init__(source_names=source_names)

    def __str__(self):
        descriptor_keys = '\n                  '.join(ii for ii in self.descriptors.keys())
        return super().__str__() + f'\n descriptor keys: {descriptor_keys}'

    def _get_uri(self, file_name, scheme):
        if self._file_uri:
            return self._file_uri
        else:
            return super()._get_uri(file_name, scheme)

    @property
    def is_fsc(self):
        return '_fsc' in self._file_name

    @property
    def is_fullframe(self):
        return '_cmos' in self._file_name

    @property
    def is_lightcurve(self):
        return self._product_id.endswith('_lcv')

    def is_valid(self):
        return True

    def set_file_id(self):
        self._file_id = TAOSIIName.remove_extensions(self._file_name)

    def set_obs_id(self):
        self._obs_id = TAOSIIName.replace_for_obs_id(self._file_id)

    @staticmethod
    def remove_extensions(value):
        return StorageName.remove_extensions(value).replace('.hdf5', '').replace('.h5', '')

    @staticmethod
    def replace_for_obs_id(value):
        return TAOSIIName.remove_extensions(value)


def set_storage_name_from_local_preconditions(storage_name, source_fqn, logger):
    """Retrieve FileInfo and header metadata into memory from files on disk. These files have extension names and
    compression as expected and support by CADC's Storage Inventory system."""
    logger.debug(f'Begin set_storage_name_from_local_preconditions for {source_fqn}')
    f_name_bits = basename(source_fqn).split('_')
    if len(f_name_bits) >= 3:
        ts = make_datetime(f_name_bits[1].replace('Z', ''))
        ts_path = f'{storage_name.scheme}:{storage_name.collection}/{ts.year}/{ts.month:02}/{ts.day:02}'
        logger.debug(f'ts_path is {ts_path}')
    else:
        raise CadcException(f'Unexpected naming pattern for file name {source_fqn}.')

    f_in = h5py.File(source_fqn)
    if storage_name.is_lightcurve:
        run_id = _lookup(f_in, '//LIGHTCURVE/HEADER/RUN/RUN_ID')
    else:
        run_id = _lookup(f_in, '//IMAGE/HEADER/RUN/RUN_ID')
    if run_id is not None:
        uri = f'{ts_path}/{run_id}/{basename(source_fqn)}'
    else:
        uri = f'{ts_path}/run_id/{basename(source_fqn)}'
    logger.info(f'Built uri {uri} from {source_fqn}')

    file_info = storage_name.file_info.get(storage_name.destination_uris[0])
    storage_name._file_uri = uri
    storage_name._destination_uris = []
    storage_name._destination_uris.append(uri)

    # now reset the storage name instance to use the destination URIs, because common code expects URIs
    storage_name._file_info = {}
    storage_name._metadata = {}
    storage_name.descriptors = {}
    storage_name._file_info[uri] = file_info
    storage_name._metadata[uri] = []
    storage_name.descriptors[uri] = f_in
    logger.debug('End set_storage_name_from_local_preconditions')


class TAOSIINoFheadVisitRunnerMeta(NoFheadVisitRunnerMeta):
    """Defines a pipeline step for all the operations that require access to the file on disk for metdata and data
    operations. This is to support HDF5 operations, since at the time of writing, there is no --fhead metadata
    retrieval option for HDF5 files.
    """

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        local_fqn = join(self._working_dir, self._storage_name.file_name)
        set_storage_name_from_local_preconditions(self._storage_name, local_fqn, self._logger)
        self._logger.debug('End _set_preconditions')


class TAOSIINoFheadLocalVisitRunnerMeta(CaomExecuteRunnerMeta):

    def __init__(self, clients, config, data_visitors, meta_visitors, reporter):
        super().__init__(clients, config, meta_visitors, reporter)
        self._data_visitors = data_visitors

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        for index, source_name in enumerate(self._storage_name.source_names):
            set_storage_name_from_local_preconditions(self._storage_name, source_name, self._logger)
        self._logger.debug('End _set_preconditions')

    def execute(self, context):
        self._logger.debug('begin execute with the steps:')
        self.storage_name = context.get('storage_name')

        self._logger.debug('set the preconditions')
        self._set_preconditions()

        self._logger.debug('get the observation for the existing model')
        self._caom2_read()

        self._logger.debug('execute the meta visitors')
        self._visit_meta()

        self._logger.debug('execute the data visitors')
        self._visit_data()

        self._logger.debug('write the observation to disk for debugging')
        self._write_model()

        self._logger.debug('store the updated xml')
        self._caom2_store()

        self._logger.debug('End execute.')


class TAOSIINoFheadScrapeVisitRunnerMeta(NoFheadScrapeRunnerMeta):

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        for source_name in self._storage_name.source_names:
            set_storage_name_from_local_preconditions(self._storage_name, source_name, self._logger)
        self._logger.debug('End _set_preconditions')


class TAOSIINoFheadStoreVisitRunnerMeta(NoFheadStoreVisitRunnerMeta):

    def __init__(self, clients, config, data_visitors, meta_visitors, reporter, store_transferrer):
        super().__init__(config, clients, store_transferrer, meta_visitors, data_visitors, reporter)

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        for source_name in self._storage_name.source_names:
            set_storage_name_from_local_preconditions(self._storage_name, source_name, self._logger)
        self._logger.debug('End _set_preconditions')


class TAOSIIOrganizeExecutesRunnerMeta(OrganizeExecutesRunnerMeta):
    """A class that extends OrganizeExecutes to handle the choosing of the correct executors based on the config.yml.
    Attributes:
        _needs_delete (bool): if True, the CAOM repo action is delete/create instead of update.
        _reporter: An instance responsible for reporting the execution status.
    Methods:
        _choose():
            Determines which descendants of CaomExecute to instantiate based on the content of the config.yml
            file for an application.
    """

    def _choose(self):
        """The logic that decides which descendants of CaomExecute to instantiate. This is based on the content of
        the config.yml file for an application.
        """
        super()._choose()
        if self._needs_delete:
            raise CadcException('No need identified for this yet.')

        if TaskType.SCRAPE in self.task_types:
            self._logger.debug(
                f'Over-riding with executor TAOSIINoFheadScrapeVisitRunnerMeta for tasks {self.task_types}.'
            )
            self._executors = []
            self._executors.append(
                TAOSIINoFheadScrapeVisitRunnerMeta(
                    self.config, self._meta_visitors, self._data_visitors, self._reporter
                )
            )
        elif TaskType.STORE in self.task_types:
            if self.config.use_local_files:
                self._logger.debug(
                    f'Over-riding with executor TAOSIINoFheadLocalVisitRunnerMeta for tasks {self.task_types}.'
                )
                self._executors = []
                self._executors.append(
                    TAOSIINoFheadStoreVisitRunnerMeta(
                        self._clients,
                        self.config,
                        self._data_visitors,
                        self._meta_visitors,
                        self._reporter,
                        self._store_transfer,
                    )
                )
            else:
                raise CadcException('Cannot store files without use_local_files set.')
        elif TaskType.MODIFY in self.task_types:
            if self.config.use_local_files:
                self._logger.debug(
                    f'Over-riding with executor TAOSIINoFheadLocalVisitRunnerMeta for tasks {self.task_types}.'
                )
                self._executors = []
                self._executors.append(
                    TAOSIINoFheadLocalVisitRunnerMeta(
                        self._clients,
                        self.config,
                        self._data_visitors,
                        self._meta_visitors,
                        self._reporter,
                    )
                )
            else:
                self._logger.debug(
                    f'Over-riding with executor TAOSIINoFheadVisitRunnerMeta for tasks {self.task_types}.'
                )
                self._executors = []
                self._executors.append(
                    TAOSIINoFheadVisitRunnerMeta(
                        self._clients,
                        self.config,
                        self._data_visitors,
                        self._meta_visitors,
                        self._modify_transfer,
                        self._reporter,
                    )
                )


class TAOSII2caom2Visitor(Fits2caom2VisitorRunnerMeta):

    def __init__(self, observation, **kwargs):
        super().__init__(observation, **kwargs)
        self._mappings = {}
        self._extension_names = defaultdict(list)

    def _get_blueprint(self, instantiated_class):
        return Hdf5ObsBlueprint(instantiated_class=instantiated_class)

    def _get_mappings(self, dest_uri):
        """
        Some file types require multiple mappings - e.g. images and lightcurves in the same file.
        """
        for key in ['IMAGE', 'LIGHTCURVE']:
            for site in ['SITE1', 'SITE2', 'SITE3', 'SITE_1', 'SITE_2', 'SITE_3']:
                # have to look for leaves - branches come back with None from _lookup
                x = f'{key}/{site}/HEADER/WCS/CNAME'
                found = _lookup(self._storage_name.descriptors.get(self._storage_name.file_uri), x)
                if found is not None:
                    self._extension_names[key].append(f'{key}/{site}')

        if self._storage_name.is_lightcurve:
            prefix = 'LIGHTCURVE'
        else:
            prefix = 'IMAGE'
        if self._storage_name.is_lightcurve:
            self._mappings[prefix] = Lightcurve(
                self._storage_name, self._clients, self._reporter, self._observation, self._config, prefix
            )
        elif '_star' in self._storage_name.file_name:
            self._mappings['IMAGE'] = TimeseriesMapping(
                self._storage_name, self._clients, self._reporter, self._observation, self._config, 'IMAGE'
            )
            if len(self._extension_names) > 1:
                # for the files that have images and lightcurves together again
                self._mappings['LIGHTCURVE'] = Lightcurve(
                    self._storage_name, self._clients, self._reporter, self._observation, self._config, 'LIGHTCURVE'
                )
        elif '_focus' in self._storage_name.file_name or '_point' in self._storage_name.file_name:
            self._mappings[prefix] = Pointing(
                self._storage_name,
                self._clients,
                self._reporter,
                self._observation,
                self._config,
                prefix,
                time_axis=3,
                energy_axis=4,
            )
        elif self._storage_name.is_fsc:
            self._mappings[prefix] = FSC(
                self._storage_name, self._clients, self._reporter, self._observation, self._config, prefix
            )
        elif (
            '_domeflat' in self._storage_name.file_name
            or '_dark' in self._storage_name.file_name
            or '_bias' in self._storage_name.file_name
        ):
            self._mappings[prefix] = DomeflatMapping(
                self._storage_name, self._clients, self._reporter, self._observation, self._config, prefix
            )
        elif self._storage_name.is_fullframe:
            self._mappings[prefix] = FullFrameImage(
                self._storage_name, self._clients, self._reporter, self._observation, self._config, prefix
            )
        else:
            self._mappings[prefix] = Single(
                self._storage_name,
                self._clients,
                self._reporter,
                self._observation,
                self._config,
                prefix,
                time_axis=3,
                energy_axis=4,
                spatial_axis_1=1,
            )

        for mapping in self._mappings.values():
            self._logger.debug(f'Created mapping {mapping.__class__.__name__}.')
        return [mapping for mapping in self._mappings.values()]

    def _get_parser(self, blueprint, uri):
        self._logger.debug(f'Using an Hdf5Parser for {self._storage_name.file_uri}')
        extension_names = []
        for key in ['IMAGE', 'LIGHTCURVE']:
            temp = blueprint._get('Observation.sequenceNumber')
            temp_lookup = None
            if isinstance(temp, tuple):
                paths, _ = temp
                temp_lookup = paths[0]
            else:
                temp_lookup = temp
            if key in temp_lookup:
                extension_names = self._extension_names.get(key)
                break

        # figure out which mapping is in use, to provide the parser with the right extension names and indices
        if extension_names is None:
            # bias, no site data
            extension_names = []
            extension_start_index = 0
            extension_end_index = 0

        self._logger.info(f'Finding data in {extension_names} for {self._storage_name.file_name}')

        if len(self._extension_names) > 1 and 'LIGHTCURVE' in extension_names[0]:
            extension_start_index=3
            extension_end_index=6  #  5 + 1 for use in range
        else:
            extension_start_index=0
            extension_end_index=3  # 2 + 1 for use in range
        # this assumes only working with one file at a time, and also, that the file is local (which, as of the time
        # of this writing, the file has to be local, until there is an --fhead option for HDF5 files)
        return Hdf5Parser(
            blueprint,
            uri,
            self._storage_name.descriptors.get(self._storage_name.file_uri),
            extension_names=extension_names,
            extension_start_index=extension_start_index,
            extension_end_index=extension_end_index,
        )


def _lookup(h5_file, things):
    """
    :params
    :things str - POSIX path breadcrumbs
    :fishing bool - True if doing a lookup for characterization, so set the logging level accordingly
    """
    things_replaced = things.replace('//', '', 1)

    def hd5f_visit(name, object):
        if isinstance(object, h5py.Dataset) and object.dtype.names is not None:
            # 'name' starts without a directory separator
            if things_replaced.startswith(name):
                for d_name in object.dtype.names:
                    temp = f'{name}/{d_name}'
                    if temp == things_replaced:
                        return object[d_name]

    result = h5_file.visititems(hd5f_visit)
    if result is not None and hasattr(result, 'decode'):
        result = result.decode('utf-8')
    return result

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

import importlib
import logging
import os
import sys
import traceback

from datetime import datetime

from astropy import wcs
from astropy.table import Table
from astropy.time import Time as AstroTime
from astropy.io.fits import Header

from caom2utils import get_gen_proc_arg_parser, gen_proc, ObsBlueprint

from caom2 import SpectralWCS, RefCoord, CoordAxis1D, Axis, CoordRange1D
from caom2 import Telescope, SimpleObservation, ObservationIntentType, Plane
from caom2 import Artifact, ProductType, ReleaseType, TemporalWCS
from caom2 import SpatialWCS, Chunk, Part, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, Dimension2D, DataProductType, Target
from caom2 import CalibrationLevel, TargetType, Provenance, Proposal
from caom2 import CoordPolygon2D, ValueCoord2D, Observation

from caom2pipe import manage_composable as mc


__all__ = ['APPLICATION', 'ARCHIVE', 'to_caom2', 'taosii_main_app',
           'COLLECTION', 'TAOSIIName']


APPLICATION = 'taosii2caom2'
COLLECTION = 'TAOSII'
ARCHIVE = 'TAOSII'


def build_energy():
    # units are nm
    min = 400
    max = 800
    central_wl = (min + max) / 2.0  # = 1200.0 / 2.0 == 600.0 nm
    fwhm = (max - min)
    ref_coord1 = RefCoord(0.5, central_wl - fwhm / 2.0)  # == 100.0 nm
    ref_coord2 = RefCoord(1.5, central_wl + fwhm / 2.0)  # == 500.0 nm
    axis = CoordAxis1D(axis=Axis(ctype='WAVE', cunit='nm'))
    axis.range = CoordRange1D(ref_coord1, ref_coord2)
    energy = SpectralWCS(axis=axis,
                         specsys='TOPOCENT',
                         ssyssrc='TOPOCENT',
                         ssysobs='TOPOCENT',
                         bandpass_name='CLEAR')
    return energy


def build_time(seconds, microseconds):
    mjd_start = AstroTime(seconds, format='unix')
    mjd_start.format = 'mjd'
    start = RefCoord(0.5, mjd_start.value)

    end_ms = seconds + (microseconds / 1e7)
    mjd_end = AstroTime(end_ms, format='unix')
    mjd_end.format = 'mjd'
    end = RefCoord(1.5, mjd_end.value)

    axis = CoordAxis1D(axis=Axis(ctype='TIME', cunit='d'))
    axis.range = CoordRange1D(start, end)
    return TemporalWCS(axis=axis,
                       timesys='UTC',
                       trefpos=None,
                       mjdref=None,
                       exposure=(microseconds / 1e7),
                       resolution=None)


def build_position(hdf5_wcs, window, telescope):
    h = Header()
    # logging.error(hdf5_wcs)
    # logging.error(window)
    h['NAXIS1'] = window['X1'].data[telescope] - window['X0'].data[telescope] + 1
    h['NAXIS2'] = window['Y1'].data[telescope] - window['Y0'].data[telescope] + 1
    h['CRVAL1'] = hdf5_wcs['CRVAL1'].data[telescope]
    h['CRVAL2'] = hdf5_wcs['CRVAL2'].data[telescope]
    h['CRPIX1'] = hdf5_wcs['CRPIX1'].data[telescope]
    h['CRPIX2'] = hdf5_wcs['CRPIX2'].data[telescope]
    # h['CDELT1'] = 4.0 / 3600.0
    # h['CDELT2'] = 4.0 / 3600.0
    #
    # JJK - use only one of CD* or CDELT* - they are either conflicting or
    # redundant
    #
    h['CD1_1'] = hdf5_wcs['CD1_1'].data[telescope]
    h['CD1_2'] = hdf5_wcs['CD1_2'].data[telescope]
    h['CD2_1'] = hdf5_wcs['CD2_1'].data[telescope]
    h['CD2_2'] = hdf5_wcs['CD2_2'].data[telescope]
    w = wcs.WCS()
    bounds = CoordPolygon2D()
    result = w.calc_footprint(header=h)
    for ii in result:
        vertex = ValueCoord2D(ii[0], ii[1])
        bounds.vertices.append(vertex)

    # ref_coord_x = RefCoord(mc.to_float(pix1), ra)
    # ref_coord_y = RefCoord(mc.to_float(pix2), dec)

    # coord = Coord2D(ref_coord_x, ref_coord_y)
    # dimension = Dimension2D(2, 2)

    # pixscale = 4.0 / 3600.0 => 4", per JJK
    # VLASS takes the cd?? approach shown here

    # function = CoordFunction2D(dimension=dimension,
    #                            ref_coord=coord,
    #                            cd11=4.0/3600.0,
    #                            cd12=0.0,
    #                            cd21=0.0,
    #                            cd22=4.0/3600.0)
    axis = CoordAxis2D(axis1=Axis(ctype='RA---TAN', cunit='deg'),
                       axis2=Axis(ctype='DEC--TAN', cunit='deg'),
                       error1=None,
                       error2=None,
                       range=None,
                       bounds=bounds,
                       function=None)

    return SpatialWCS(axis=axis,
                      coordsys='FK5',
                      equinox=2000.0,
                      resolution=None)


def build_from_hdf5(args):
    obs = None
    if args.in_obs_xml:
        obs = mc.read_obs_from_file(args.in_obs_xml.name)
    index = 0
    for f_name in args.local:
        product_id = args.lineage[index].split('/')[0]
        t_header = Table.read(f_name, format='hdf5', path='header')
        # logging.error(t_header.colnames)
        # ['VERSION_MAJOR', 'VERSION_MINOR', 'TIME_IN_SEC',
        # 'TIME_IN_MICROSEC', 'RUN_ID', 'ORIGIN', 'OBSMODE', 'FIELD', 'RA',
        # 'DEC', 'EXPTIME', 'NUM_IMAGER']
        # logging.error(t_header)

        # t_header['RUN_ID'].data[0].decode() - return this string
        # 20190805T024026
        # logging.error(t_header['RUN_ID'].data[0].decode())
        release_date = datetime.strptime(
            t_header['RUN_ID'].data[0].decode(), '%Y%m%dT%H%M%S')

        t_image = Table.read(f_name, format='hdf5', path='image')
        # logging.error(t_image.colnames)
        # ['col0', 'col1', 'col2']
        # logging.error(t_image)
        # logging.error(t_image['col0'].data[0])
        # logging.error(t_image['col0'].data[143999])

        t_catalog = Table.read(f_name, format='hdf5', path='catalog')
        # logging.error(t_catalog.colnames)
        # ['CAT_ID', 'GAIA_ID', '2MASS_ID', 'RA', 'DEC', 'TAOS_MAG',
        # 'GAIA_MAG', '2MASS_JMAG']
        # logging.error(t_catalog)

        t_imager = Table.read(f_name, format='hdf5', path='imager')
        # logging.error(t_imager.colnames)
        # ['TEL_ID', 'CAM_ID', 'IMGR_ID', 'XLOC', 'YLOC']
        # logging.error(t_imager)

        t_moment = Table.read(f_name, format='hdf5', path='moment')
        # logging.error(t_moment.colnames)
        # ['col0', 'col1', 'col2']
        # logging.error(t_moment)

        t_window = Table.read(f_name, format='hdf5', path='window')
        # logging.error(t_window.colnames)
        # ['X0', 'X1', 'Y0', 'Y1', 'XC', 'YC']
        # logging.error(t_window)

        t_wcs= Table.read(f_name, format='hdf5', path='/wcs/cdmatrix')
        # ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2']
        # logging.error(t_wcs)

        taos = Telescope(name='TAOS',
                         geo_location_x=-2354953.99637757,
                         geo_location_y=-4940160.3636381,
                         geo_location_z=3270123.70695983)

        target = Target(name=str(t_header['FIELD'].data[0]),
                        target_type=TargetType.FIELD,
                        standard=None,
                        redshift=None,
                        keywords=None,
                        moving=None)

        proposal = Proposal(id=COLLECTION,
                            pi_name=None,
                            project=COLLECTION,
                            title=None)

        if obs is None:
            obs = SimpleObservation(collection=COLLECTION,
                                    observation_id=args.observation[1],
                                    sequence_number=None,
                                    intent=ObservationIntentType.SCIENCE,
                                    type='FIELD',
                                    proposal=proposal,
                                    telescope=taos,
                                    instrument=None,
                                    target=target,
                                    meta_release=release_date)

        provenance = Provenance(name=COLLECTION,
                                version='{}.{}'.format(
                                    t_header['VERSION_MAJOR'].data[0],
                                    t_header['VERSION_MINOR'].data[0]),
                                project=COLLECTION,
                                producer=COLLECTION,
                                run_id=t_header['RUN_ID'].data[0].decode(),
                                reference='https://taos2.asiaa.sinica.edu.tw/',
                                last_executed=release_date)

        plane = Plane(product_id=product_id,
                      data_release=release_date,
                      meta_release=release_date,
                      provenance=provenance,
                      data_product_type=DataProductType.IMAGE,
                      calibration_level=CalibrationLevel.RAW_STANDARD)

        artifact = mc.get_artifact_metadata(
            f_name, ProductType.SCIENCE, ReleaseType.DATA,
            mc.build_uri(COLLECTION, os.path.basename(f_name)))

        # parts are always named '0'
        part = Part('0')

        # do each of the three telescopes
        for telescope in [0, 1, 2]:
            position = build_position(t_wcs,
                                      t_window,
                                      telescope)

            time = build_time(t_header['TIME_IN_SEC'].data[0],
                              t_header['TIME_IN_MICROSEC'].data[0])

            energy = build_energy()

            chunk = Chunk(naxis=4,
                          position_axis_1=1,
                          position_axis_2=2,
                          energy_axis=3,
                          time_axis=4,
                          position=position,
                          energy=energy,
                          time=time)

            part.chunks.append(chunk)

        artifact.parts.add(part)
        plane.artifacts.add(artifact)
        obs.planes.add(plane)

        index += 1

    return obs


# bash-5.0#  /usr/local/hdf5/bin/h5dump 20190805T024026_f060_s00001.h5 | grep DATASET
#    DATASET "catalog" {
#    DATASET "header" {
#    DATASET "image" {
#    DATASET "imager" {
#    DATASET "moment" {
#       DATASET "a" {
#       DATASET "ap" {
#       DATASET "b" {
#       DATASET "bp" {
#       DATASET "cdmatrix" {
#       DATASET "order" {
#    DATASET "window" {


class TAOSIIName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage
    """

    TAOSII_NAME_PATTERN = '*'

    def __init__(self, file_name=None, entry=None):
        self._file_name = file_name
        self._obs_id = TAOSIIName.get_obs_id(file_name)
        super(TAOSIIName, self).__init__(
            self._obs_id, COLLECTION, TAOSIIName.TAOSII_NAME_PATTERN,
            file_name, entry=entry)

    @property
    def file_name(self):
        return self._file_name

    @property
    def product_id(self):
        return 'pixel'

    def is_valid(self):
        return True

    @staticmethod
    def get_obs_id(value):
        bits = value.split('.')[0].split('_')
        return f'{bits[1]}_{bits[2]}'


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    bp.configure_position_axes((1,2))
    bp.configure_time_axis(3)
    bp.configure_energy_axis(4)
    bp.configure_polarization_axis(5)
    bp.configure_observable_axis(6)
    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes (an n:n
    relationship between TDM attributes and CAOM attributes). Must have this
    signature for import_module loading and execution.
    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = kwargs.get('headers')
    fqn = kwargs.get('fqn')
    uri = kwargs.get('uri')
    taosii_name = None
    if uri is not None:
        taosii_name = TAOSIIName(artifact_uri=uri)
    if fqn is not None:
        taosii_name = TAOSIIName(file_name=os.path.basename(fqn))
    if taosii_name is None:
        raise mc.CadcException(f'Need one of fqn or uri defined for '
                               f'{observation.observation_id}')
    logging.debug('Done update.')
    return observation


def _build_blueprints(uris):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.
    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.
    :param uris The artifact URIs for the files to be processed."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not mc.StorageName.is_preview(uri):
            accumulate_bp(blueprint, uri)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.local:
        for ii in args.local:
            file_id = mc.StorageName.remove_extensions(os.path.basename(ii))
            file_name = f'{file_id}.fits'
            result.append(TAOSIIName(file_name=file_name).file_uri)
    elif args.lineage:
        for ii in args.lineage:
            result.append(ii.split('/', 1)[1])
    else:
        raise mc.CadcException(
            f'Could not define uri from these args {args}')
    return result


def to_caom2():
    """This function is called by pipeline execution. It must have this name.
    """
    args = get_gen_proc_arg_parser().parse_args()
    # uris = _get_uris(args)
    # blueprints = _build_blueprints(uris)
    # result = gen_proc(args, blueprints)
    obs = build_from_hdf5(args)
    mc.write_obs_to_file(obs, f'./{obs.observation_id}.actual.xml')
    logging.debug(f'Done {APPLICATION} processing.')
    return 1


def taosii_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = to_caom2()
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)

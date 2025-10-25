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

import logging
from taosii2caom2 import set_storage_name_from_local_preconditions, TAOSIIName


def test_is_valid():
    assert TAOSIIName(['20190805T024026_f060_s00001']).is_valid()


def test_storage_name(test_config, test_data_dir):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    for test_f_id in [
        'taos2_20240221T014139Z_star987654321',
        'taos2_20240221T033316Z_cmos31_012',
        'taos2_20240221T033325Z_star987654321',
        'taos2_20240221T033329Z_star987654321_lcv',
        'taos2_20241116T092303Z_cmos31_bias_012',
        'taos2_20241116T092310Z_cmos31_012',
        'taos2_20241116T092321Z_star987654321',
        'taos2_20241121T100428Z_fsc_145',
        'taos2_20241121T100920Z_cmos31_012',
        'taos2_20241121T101011Z_star987654321',
        'taos2_20250429T032413Z_cmos00_skyflat_001',
    ]:
        for prefix in ['', '/data/', 'cadc:TAOSII/']:
            test_f_name = f'{test_f_id}.h5'
            source_names = [f'{prefix}{test_f_id}.h5']
            test_subject = TAOSIIName(source_names=source_names)
            source_fqn = f'{test_data_dir}/2024/try2/{test_f_name}'
            set_storage_name_from_local_preconditions(test_subject, source_fqn, logger)
            assert test_subject.obs_id == test_f_id, f'wrong obs id {test_subject.obs_id}'
            assert test_subject.file_name == test_f_name, 'wrong file name'
            assert test_subject.product_id == test_f_id, 'wrong product id'
            assert len(test_subject.destination_uris) == 1, 'destination uris not set'
            if '_lcv' in test_f_id:
                assert test_subject.is_lightcurve, f'lightcurve {test_f_id}'
            else:
                assert not test_subject.is_lightcurve, f'not lightcurve {test_f_id}'
            if not (
                test_subject.file_uri.startswith('cadc:TAOSII/2024/02/21/')
                or test_subject.file_uri.startswith('cadc:TAOSII/2024/11/16/')
                or test_subject.file_uri.startswith('cadc:TAOSII/2024/11/21/')
                or test_subject.file_uri.startswith('cadc:TAOSII/2025/04/29/')
            ):
                assert False, f'wrong file_uri {test_subject.file_uri}'

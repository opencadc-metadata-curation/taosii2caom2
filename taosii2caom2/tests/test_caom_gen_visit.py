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
import traceback

from os import unlink
from os.path import dirname, exists, join, realpath

from caom2utils.data_util import get_local_file_info
from caom2.diff import get_differences
from caom2pipe.manage_composable import ExecutionReporter2, read_obs_from_file, write_obs_to_file
from taosii2caom2 import set_storage_name_from_local_preconditions, TAOSIIName
from taosii2caom2 import file2caom2_augmentation

from glob import glob

THIS_DIR = dirname(realpath(__file__))
TEST_DATA_DIR = join(THIS_DIR, 'data')
PLUGIN = join(dirname(THIS_DIR), 'main_app.py')


def pytest_generate_tests(metafunc):
    obs_id_list = glob(f'{TEST_DATA_DIR}/2024/try2/*.h5')
    metafunc.parametrize('test_name', obs_id_list)


def test_visitor(test_config, test_name, tmp_path):
    """
    taos2_20240208T233041Z_cmos31_012.h5 is a full-frame image file
    taos2_20240208T233045Z_star987654321.h5 is a window mode image file
    taos2_20240208T233034Z_star987654321_lcv.h5 is a lightcurve file
    taos2_20240208T232951Z_star987654321.h5 is a combined window mode image/lightcurve file
    taos2_20240209T191043Z_fsc_145.h5 is a finder scope image file
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    test_config.change_working_directory(tmp_path.as_posix())
    storage_name = TAOSIIName(source_names=[test_name])
    storage_name.file_info[storage_name.file_uri] = get_local_file_info(test_name)
    set_storage_name_from_local_preconditions(storage_name, test_name, logger)
    test_reporter = ExecutionReporter2(test_config)
    kwargs = {
        'storage_name': storage_name,
        'config': test_config,
        'reporter': test_reporter,
    }
    expected_fqn = f'{dirname(test_name)}/{storage_name.file_id}.expected.xml'
    actual_fqn = expected_fqn.replace('expected', 'actual')
    if exists(actual_fqn):
        unlink(actual_fqn)
    observation = None
    try:
        observation = file2caom2_augmentation.visit(observation, **kwargs)
    except Exception as e:
        logging.error(e)
        logging.error(traceback.format_exc())
        assert False, f'{e}'

    if exists(expected_fqn):
        expected = read_obs_from_file(expected_fqn)
        compare_result = get_differences(expected, observation)
        if compare_result is not None:
            if observation is None:
                assert False, f'observation should not be None {test_name}'
            else:
                write_obs_to_file(observation, actual_fqn)
                compare_text = '\n'.join([r for r in compare_result])
                msg = f'Differences found in observation {expected.observation_id}\n{compare_text}'
                raise AssertionError(msg)
    else:
        if observation:
            write_obs_to_file(observation, actual_fqn)
            assert False, 'no expected observation'
        else:
            assert False, 'no Observation to check against, and no expected Observation to check with.'
    # assert False


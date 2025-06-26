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

from collections import deque
from datetime import datetime, timedelta
from mock import patch, PropertyMock

from caom2pipe.data_source_composable import RunnerMeta
from caom2pipe.manage_composable import Config, State, TaskType
from taosii2caom2 import composable, TAOSIIName

TEST_TIME = datetime(year=2019, month=3, day=7, hour=19, minute=5)


@patch('caom2pipe.data_source_composable.QueryTimeBoxDataSourceRunnerMeta._initialize_end_dt')
@patch(
    'caom2pipe.data_source_composable.QueryTimeBoxDataSourceRunnerMeta.end_dt',
    new_callable=PropertyMock(return_value=datetime(year=2019, month=3, day=7, hour=19, minute=5))
)
@patch('caom2pipe.client_composable.ClientCollection')
@patch(
    'caom2pipe.data_source_composable.QueryTimeBoxDataSourceRunnerMeta.get_time_box_work', autospec=True,
)
@patch('caom2pipe.execute_composable.OrganizeExecutesRunnerMeta.do_one')
def test_run_incremental(
    run_mock,
    get_work_mock,
    clients_mock,
    end_time_mock,
    initialize_end_dt_mock,
    test_data_dir,
    test_config,
    tmp_path,
    change_test_dir,
):
    test_obs_id = 'taos2_20240221T033329Z_star987654321_lcv'
    test_f_name = f'{test_obs_id}.h5'
    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.interval = 7200
    test_config.task_types = [TaskType.MODIFY]
    test_config.logging_level = 'DEBUG'

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    test_state_fqn = f'{tmp_path}/state.yml'
    start_time = datetime(year=2019, month=3, day=3, hour=19, minute=5)
    State.write_bookmark(test_state_fqn, test_config.bookmark, start_time)
    Config.write_to_file(test_config)
    run_mock.return_value = (0, None)

    def _mock_query_listing(_, ignore_prev_exec_dt, ignore_exec_dt):
        result = deque()
        dt = TEST_TIME - timedelta(minutes=5)
        storage_name = TAOSIIName(
            source_names=[f'{test_data_dir}/2024/try2/taos2_20240221T033329Z_star987654321_lcv.h5']
        )
        result.append(RunnerMeta(storage_name, dt))
        return result

    get_work_mock.side_effect = _mock_query_listing

    # execution
    test_result = composable._run_incremental()
    assert test_result == 0, 'mocking correct execution'
    assert run_mock.called, 'should have been called'
    args, _ = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, TAOSIIName), type(test_storage)
    assert test_storage.obs_id == test_obs_id, f'wrong obs id {test_storage.obs_id}'
    assert test_storage.file_name == test_f_name, 'wrong file name'
    assert test_storage.file_uri == f'cadc:TAOSII/{test_f_name}', 'wrong uri'
    assert not clients_mock.data_client.get.called, 'data get'
    assert not clients_mock.data_client.get_head.called, 'data get_head'
    assert not clients_mock.data_client.info.called, 'data info'
    assert not clients_mock.data_client.put.called, 'data put'
    assert not clients_mock.data_client.remove.called, 'data remove'
    assert not clients_mock.metadata_client.create.called, 'metadata create'
    assert not clients_mock.metadata_client.read.called, 'metadata read'
    assert not clients_mock.metadata_client.update.called, 'metadata update'
    assert not clients_mock.metadata_client.delete.called, 'metadata delete'


@patch('caom2pipe.client_composable.ClientCollection')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run(run_mock, clients_mock, test_config, tmp_path, change_test_dir):
    test_f_id = 'test_file_id'
    test_f_name = f'{test_f_id}.hdf5'
    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.MODIFY]

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    with open(test_config.work_fqn, 'w') as f:
        f.write(test_f_name)

    Config.write_to_file(test_config)

    composable._run()
    assert run_mock.called, 'should have been called'
    args, _ = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, TAOSIIName), type(test_storage)
    assert test_storage.obs_id == test_f_id, 'wrong obs id'
    assert test_storage.file_name == test_f_name, 'wrong file name'


@patch('caom2pipe.client_composable.ClientCollection')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run_store_modify(run_mock, clients_mock, test_config, tmp_path, change_test_dir):
    test_config.change_working_directory(tmp_path)
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.STORE, TaskType.MODIFY]
    test_config.use_local_files = True
    test_config.data_sources = ['/usr/src/app/taosii2caom2/int_test/test_files']
    test_config.data_source_extensions = ['.h5']
    test_config.store_modified_files_only = True
    test_config.logging_level = 'DEBUG'
    run_mock.return_value = (0, None)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    Config.write_to_file(test_config)

    composable._run()
    assert run_mock.called, 'should have been called'
    assert run_mock.call_count == 3, 'wrong number of calls'
    args, _ = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, TAOSIIName), type(test_storage)
    assert test_storage.obs_id == 'taos2_20241116T092321Z_star987654321', f'wrong obs id {test_storage}'

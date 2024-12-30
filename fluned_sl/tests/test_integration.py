"""
fluned_sl testing
"""
import os
from fluned_sl.fluned_sl_case import FlunedSlCase

TEST_REL_PATH_CASE_1 = "test/case_1"
TEST_REL_PATH_CASE_2 = "test/case_2"
TEST_REL_PATH_CASE_3 = "test/case_3"
TEST_REL_PATH_CASE_4 = "test/case_4"
TEST_REL_PATH_CASE_5 = "test/case_5"
TEST_REL_PATH_CASE_6 = "test/case_6"
TEST_REL_PATH_CASE_7 = "test/case_7"
TEST_REL_PATH_CASE_8 = "test/case_8"
TEST_REL_PATH_CASE_9 = "test/case_9"


def isclose(a, b, rel_tol=1e-04, abs_tol=1e-08):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def test_tot_inflows_case_2():
    """
    This function tests the total inflows
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_2)
    test_case = FlunedSlCase(test_setting, test_path)
    assert len(test_case.get_inflow_node_list()) == 1

def test_tot_inflow_mass_case_2():
    """
    This function tests the total inflowing mass
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_2)
    test_case = FlunedSlCase(test_setting, test_path)
    assert test_case.get_total_inflow() == 10

def test_tot_outflows_case_2():
    """
    This function check the total number of outflow nodes
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_2)
    test_case = FlunedSlCase(test_setting, test_path)
    assert len(test_case.get_outflow_node_list()) == 2

def test_tot_outflow_mass_case_2():
    """
    This function tests the total outflowing mass
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_2)
    test_case = FlunedSlCase(test_setting, test_path)
    assert test_case.get_total_outflow() == 10

def test_inflow_outflow_balance_case_2():
    """
    This function check that the inflow and outflow are balanced
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_2)
    test_case = FlunedSlCase(test_setting, test_path)
    assert test_case.get_total_outflow() == test_case.get_total_inflow()

def test_inflow_outflow_balance_case_3():
    """
    This function check that the inflow and outflow are balanced
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    assert test_case.get_total_outflow() == test_case.get_total_inflow()

def test_density_case_3():
    """
    This function check the density values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('density_g_cm3')
    expected_values = {('01_DECAY_TEST', 1):0.99426,
                       ('01_DECAY_TEST', 2):0.99426,
                       ('01_DECAY_TEST', 3):0.99426,
                       ('01_DECAY_TEST', 4):0.99426,
                       ('01_DECAY_TEST', 5):0.99426,
                       }
    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_reynolds_case_3():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('reynolds')
    expected_values = {('01_DECAY_TEST', 1):91796.76895,
                       ('01_DECAY_TEST', 2):30598.92298,
                       ('01_DECAY_TEST', 3):7649.730746,
                       ('01_DECAY_TEST', 4):917.9676895,
                       ('01_DECAY_TEST', 5):18359.35379,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=0, atol=1)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_bulk_velocity_case_3():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('bulk_velocity_m_s')
    expected_values = {('01_DECAY_TEST', 1):3.201475330,
                       ('01_DECAY_TEST', 2):0.355719481,
                       ('01_DECAY_TEST', 3):0.022232468,
                       ('01_DECAY_TEST', 4):0.000320148,
                       ('01_DECAY_TEST', 5):0.128059013,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-03, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_bulk_res_time_case_3():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('bulk_res_time')
    expected_values = {('01_DECAY_TEST', 1):0.312355991,
                       ('01_DECAY_TEST', 2):2.811203921,
                       ('01_DECAY_TEST', 3):44.97926273,
                       ('01_DECAY_TEST', 4):3123.559912,
                       ('01_DECAY_TEST', 5):23.42669934,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-03, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_res_time_pipe_uniform_case_3():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_res_time')
    expected_values = {('01_DECAY_TEST', 1):0.312355991,
                       ('01_DECAY_TEST', 2):2.811203921,
                       ('01_DECAY_TEST', 3):44.97926273,
                       ('01_DECAY_TEST', 4):3123.559912,
                       ('01_DECAY_TEST', 5):23.42669934,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_3():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-03
    test_setting = {
            "pipe_time_default": "uniform",
             "mc_error_max" : precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):0.931299741,
                       ('01_DECAY_TEST', 4):1.2353E-132,
                       ('01_DECAY_TEST', 5):1.2668E-133,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_3_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-03
    test_setting = {
            "pipe_time_default": "uniform",
             "mc_error_max" : precision,
             "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):0.931299741,
                       ('01_DECAY_TEST', 4):1.2353E-132,
                       ('01_DECAY_TEST', 5):1.2668E-133,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)


def test_av_activity_pipe_uniform_case_3():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):98.49695854,
                       ('01_DECAY_TEST', 2):84.88250417,
                       ('01_DECAY_TEST', 3):16.66709666,
                       ('01_DECAY_TEST', 4):0.003066929,
                       ('01_DECAY_TEST', 5):4.868E-133 ,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_av_activity_pipe_uniform_case_3_det():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
            "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):98.49695854,
                       ('01_DECAY_TEST', 2):84.88250417,
                       ('01_DECAY_TEST', 3):16.66709666,
                       ('01_DECAY_TEST', 4):0.003066929,
                       ('01_DECAY_TEST', 5):4.868E-133 ,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_tot_activity_pipe_uniform_case_3():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('tot_activity_bq')
    expected_values = {('01_DECAY_TEST', 1):0.030943732,
                       ('01_DECAY_TEST', 2):0.239999626,
                       ('01_DECAY_TEST', 3):0.754001689,
                       ('01_DECAY_TEST', 4):0.009635042,
                       ('01_DECAY_TEST', 5):1.147E-134 ,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_tot_activity_pipe_uniform_case_3_det():
    """
    This function check the reynolds values against calculated ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
            "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('tot_activity_bq')
    expected_values = {('01_DECAY_TEST', 1):0.030943732,
                       ('01_DECAY_TEST', 2):0.239999626,
                       ('01_DECAY_TEST', 3):0.754001689,
                       ('01_DECAY_TEST', 4):0.009635042,
                       ('01_DECAY_TEST', 5):1.147E-134 ,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_cumulated_res_time_pipe_uniform_case_3():
    """
    This function check the cumulated residence time values against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_3)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_time_cumulated')
    expected_values = {('01_DECAY_TEST', 1):0.312355991,
                       ('01_DECAY_TEST', 2):3.123559912,
                       ('01_DECAY_TEST', 3):48.10282264,
                       ('01_DECAY_TEST', 4):3171.662734,
                       ('01_DECAY_TEST', 5):3195.089434,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_inflow_outflow_balance_case_4():
    """
    This function check that the inflow and outflow are balanced
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    assert test_case.get_total_outflow() == test_case.get_total_inflow()

def test_density_case_4():
    """
    This function check the density values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('density_g_cm3')
    expected_values = {('branches', 1):0.95743,
                       ('branches', 2):0.95743,
                       ('branches', 3):0.95743,
                       ('branches', 4):0.98561,
                       ('branches', 5):0.98561,
                       ('branches', 6):0.95743,
                       ('branches', 7):0.95743,
                       ('branches', 8):0.98561,
                       ('branches', 9):0.98561,
                       ('branches',10):0.95743,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_mass_flow_case_4():
    """
    This function check the mass flow values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('mass_flow')
    expected_values = {('branches', 1):20,
                       ('branches', 2):20,
                       ('branches', 3):3 ,
                       ('branches', 4):14,
                       ('branches', 5):14,
                       ('branches', 6):1 ,
                       ('branches', 7):2 ,
                       ('branches', 8):10,
                       ('branches', 9):18,
                       ('branches',10):12,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-05, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)


def test_bulk_velocity_case_4():
    """
    This function check the bulk velocity values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('bulk_velocity_m_s')
    expected_values = {('branches', 1):4.380616101,
                       ('branches', 2):4.380616101,
                       ('branches', 3):0.657092415,
                       ('branches', 4):2.978757613,
                       ('branches', 5):2.978757613,
                       ('branches', 6):0.219030805,
                       ('branches', 7):0.43806161 ,
                       ('branches', 8):2.12768401 ,
                       ('branches', 9):3.829831217,
                       ('branches',10):2.628369661,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_bulk_res_time_case_4():
    """
    This function check the bulk residence time values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('bulk_res_time')
    expected_values = {('branches', 1):0.164585525,
                       ('branches', 2):0.164585525,
                       ('branches', 3):2.194470621,
                       ('branches', 4):0.484084369,
                       ('branches', 5):0.484084369,
                       ('branches', 6):6.583411862,
                       ('branches', 7):3.291705931,
                       ('branches', 8):1.382710020,
                       ('branches', 9):0.376510065,
                       ('branches',10):0.548617655,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check = isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_bulk_reynolds_case_4():
    """
    This function check the bulk reynolds values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('reynolds')
    expected_values = {('branches', 1):1182155.42,
                       ('branches', 2):1182155.42,
                       ('branches', 3):177323.31 ,
                       ('branches', 4):467257.92 ,
                       ('branches', 5):467257.92 ,
                       ('branches', 6):59107.77  ,
                       ('branches', 7):118215.54 ,
                       ('branches', 8):333755.66 ,
                       ('branches', 9):600760.18 ,
                       ('branches',10):709293.25 ,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_4():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('branches', 1):98.41270456,
                       ('branches', 2):96.85060418,
                       ('branches', 3):78.24405656,
                       ('branches', 4):92.39834624,
                       ('branches', 5):88.15076023,
                       ('branches', 6):51.06806546,
                       ('branches', 7):70.32761189,
                       ('branches', 8):524.5331793,
                       ('branches', 9):81.40466439,
                       ('branches',10):425.5213079,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_4_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
            "numerical_method" : "deterministic"
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('branches', 1):98.41270456,
                       ('branches', 2):96.85060418,
                       ('branches', 3):78.24405656,
                       ('branches', 4):92.39834624,
                       ('branches', 5):88.15076023,
                       ('branches', 6):51.06806546,
                       ('branches', 7):70.32761189,
                       ('branches', 8):524.5331793,
                       ('branches', 9):81.40466439,
                       ('branches',10):425.5213079,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_cumulated_res_time_uniform_case_4():
    """
    this function check the cumulated residence time in each pipe against
    calculated ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_4)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_time_cumulated')

    expected_values = {('branches', 1):0.164585525,
                       ('branches', 2):0.329171050,
                       ('branches', 3):2.523641670,
                       ('branches', 4):0.813255419,
                       ('branches', 5):1.297339788,
                       ('branches', 6):6.912582911,
                       ('branches', 7):3.620876981,
                       ('branches', 8):1.382710020,
                       ('branches', 9):2.190191451,
                       ('branches',10):2.304355502,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)


def test_mass_flow_case_5():
    """
    This function check the mass flow values against calculated ones
    """
    test_setting = {}
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_dict = test_case.get_all_attr('mass_flow')
    expected_values = {('01_CIRCUIT_A', 1):20,
                       ('01_CIRCUIT_A', 2):20,
                       ('03_CIRCUIT_C', 1):3 ,
                       ('03_CIRCUIT_C', 2):14 ,
                       ('03_CIRCUIT_C', 3):18 ,
                       ('04_CIRCUIT_D', 1):14,
                       ('02_CIRCUIT_B', 1):1 ,
                       ('02_CIRCUIT_B', 2):2 ,
                       ('05_CIRCUIT_E', 1):10 ,
                       ('05_CIRCUIT_E', 2):12 ,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-05, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_5():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('01_CIRCUIT_A', 1):98.41270456,
                       ('01_CIRCUIT_A', 2):96.85060418,
                       ('03_CIRCUIT_C', 1):78.24405656,
                       ('03_CIRCUIT_C', 2):88.15076023,
                       ('03_CIRCUIT_C', 3):81.40466439,
                       ('04_CIRCUIT_D', 1):92.39834624,
                       ('02_CIRCUIT_B', 1):51.06806546,
                       ('02_CIRCUIT_B', 2):70.32761189,
                       ('05_CIRCUIT_E', 1):524.5331793,
                       ('05_CIRCUIT_E', 2):425.5213079,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_out_activity_pipe_uniform_case_5_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
             "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('01_CIRCUIT_A', 1):98.41270456,
                       ('01_CIRCUIT_A', 2):96.85060418,
                       ('03_CIRCUIT_C', 1):78.24405656,
                       ('03_CIRCUIT_C', 2):88.15076023,
                       ('03_CIRCUIT_C', 3):81.40466439,
                       ('04_CIRCUIT_D', 1):92.39834624,
                       ('02_CIRCUIT_B', 1):51.06806546,
                       ('02_CIRCUIT_B', 2):70.32761189,
                       ('05_CIRCUIT_E', 1):524.5331793,
                       ('05_CIRCUIT_E', 2):425.5213079,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_out_activity_pipe_uniform_transient_10s_case_5():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
            "t_bins":"10"
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()

    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('01_CIRCUIT_A', 1):98.41270456,
                       ('01_CIRCUIT_A', 2):96.85060418,
                       ('03_CIRCUIT_C', 1):78.24405656,
                       ('03_CIRCUIT_C', 2):88.15076023,
                       ('03_CIRCUIT_C', 3):81.40466439,
                       ('04_CIRCUIT_D', 1):92.39834624,
                       ('02_CIRCUIT_B', 1):51.06806546,
                       ('02_CIRCUIT_B', 2):70.32761189,
                       ('05_CIRCUIT_E', 1):524.5331793,
                       ('05_CIRCUIT_E', 2):425.5213079,
                       }


    a_vec = []
    b_vec = []
    check_vec = []


    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_out_activity_pipe_uniform_transient_01s_case_5():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
            "t_bins":"0.1"
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('01_CIRCUIT_A', 1):0,
                       ('01_CIRCUIT_A', 2):0,
                       ('03_CIRCUIT_C', 1):0,
                       ('03_CIRCUIT_C', 2):0,
                       ('03_CIRCUIT_C', 3):0,
                       ('04_CIRCUIT_D', 1):0,
                       ('02_CIRCUIT_B', 1):0,
                       ('02_CIRCUIT_B', 2):0,
                       ('05_CIRCUIT_E', 1):0,
                       ('05_CIRCUIT_E', 2):0,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_average_activity_pipe_uniform_transient_01s_case_5():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-3
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
            "t_bins":"0.1"
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')

    expected_values = {('01_CIRCUIT_A', 1):60.27518882,
                       ('01_CIRCUIT_A', 2):0,
                       ('03_CIRCUIT_C', 1):0,
                       ('03_CIRCUIT_C', 2):0,
                       ('03_CIRCUIT_C', 3):0,
                       ('04_CIRCUIT_D', 1):0,
                       ('02_CIRCUIT_B', 1):0,
                       ('02_CIRCUIT_B', 2):0,
                       ('05_CIRCUIT_E', 1):43.18280417,
                       ('05_CIRCUIT_E', 2):0,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_average_activity_pipe_uniform_transient_1s_case_5():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-3
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
            "t_bins":"1"
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_5)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')

    expected_values = {('01_CIRCUIT_A', 1):99.20423586,
                       ('01_CIRCUIT_A', 2):97.62957154,
                       ('03_CIRCUIT_C', 1):28.66157512,
                       ('03_CIRCUIT_C', 2):35.32278315,
                       ('03_CIRCUIT_C', 3):0,
                       ('04_CIRCUIT_D', 1):94.60701538,
                       ('02_CIRCUIT_B', 1):9.553858375,
                       ('02_CIRCUIT_B', 2):19.10771675,
                       ('05_CIRCUIT_E', 1):413.50527,
                       ('05_CIRCUIT_E', 2):0,
                       }




    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_out_activity_pipe_uniform_case_6():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_6)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):49.23639931,
                       ('01_DECAY_TEST', 4):23.04586531,
                       ('01_DECAY_TEST', 5):2.363269345,
                       ('01_DECAY_TEST', 6):1.106165109,
                       ('01_DECAY_TEST', 7):0.517757847,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_6_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
             "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_6)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):49.23639931,
                       ('01_DECAY_TEST', 4):23.04586531,
                       ('01_DECAY_TEST', 5):2.363269345,
                       ('01_DECAY_TEST', 6):1.106165109,
                       ('01_DECAY_TEST', 7):0.517757847,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_7():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):61.41676733,
                       ('01_DECAY_TEST', 4):28.7470767 ,
                       ('01_DECAY_TEST', 5):2.947907758,
                       ('01_DECAY_TEST', 6):1.379814245,
                       ('01_DECAY_TEST', 7):0.645843597,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_7_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
             "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):61.41676733,
                       ('01_DECAY_TEST', 4):28.7470767 ,
                       ('01_DECAY_TEST', 5):2.947907758,
                       ('01_DECAY_TEST', 6):1.379814245,
                       ('01_DECAY_TEST', 7):0.645843597,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)


def test_average_activity_pipe_uniform_case_7():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):98.49695854,
                       ('01_DECAY_TEST', 2):84.88250417,
                       ('01_DECAY_TEST', 3):67.41517448,
                       ('01_DECAY_TEST', 4):43.03474725,
                       ('01_DECAY_TEST', 5):11.32814232,
                       ('01_DECAY_TEST', 6):2.065599848,
                       ('01_DECAY_TEST', 7):0.966836254,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_average_activity_pipe_uniform_case_7_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    test_setting = {
            "pipe_time_default": "uniform",
            "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):98.49695854,
                       ('01_DECAY_TEST', 2):84.88250417,
                       ('01_DECAY_TEST', 3):67.41517448,
                       ('01_DECAY_TEST', 4):43.03474725,
                       ('01_DECAY_TEST', 5):11.32814232,
                       ('01_DECAY_TEST', 6):2.065599848,
                       ('01_DECAY_TEST', 7):0.966836254,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=1e-04, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_7_100s():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "t_bins":"100",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):61.41676733,
                       ('01_DECAY_TEST', 4):28.7470767 ,
                       ('01_DECAY_TEST', 5):2.947907758,
                       ('01_DECAY_TEST', 6):1.379814245,
                       ('01_DECAY_TEST', 7):0.645843597,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_7_5s():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "t_bins":"5",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):36.3392683,
                       ('01_DECAY_TEST', 4):0,
                       ('01_DECAY_TEST', 5):0,
                       ('01_DECAY_TEST', 6):0,
                       ('01_DECAY_TEST', 7):0,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_average_activity_pipe_uniform_case_7_5s():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "t_bins":"5",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
    # the values for nodes 3 and 4 have been simulated, as it is not possible
    # to calculate them analytically. The important is that the rtd
    # function correctly by passing a small fraction of the activity to the
    # following node for every time step
    expected_values = {('01_DECAY_TEST', 1):98.49695854,
                       ('01_DECAY_TEST', 2):84.88250417,
                       ('01_DECAY_TEST', 3):63.56114903,
                       ('01_DECAY_TEST', 4):1.9959910575667026,
                       ('01_DECAY_TEST', 5):0,
                       ('01_DECAY_TEST', 6):0,
                       ('01_DECAY_TEST', 7):0,
                       }




    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_7_10s():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "t_bins":"10",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):97.00905397,
                       ('01_DECAY_TEST', 2):73.81127268,
                       ('01_DECAY_TEST', 3):61.416703,
                       ('01_DECAY_TEST', 4):0.0,
                       ('01_DECAY_TEST', 5):0.0,
                       ('01_DECAY_TEST', 6):0.0,
                       ('01_DECAY_TEST', 7):0.0,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_average_activity_pipe_uniform_case_7_10s():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "steady_state": False,
            "pipe_time_default": "uniform",
            "t_bins":"10",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_7)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')

    # the values for nodes 3 and 4 have been simulated, as it is not possible
    # to calculate them analytically. The important is that the rtd
    # function correctly passes a small fraction of the activity to the
    # following node for every time step

    expected_values = {('01_DECAY_TEST', 1):98.49695854,
                       ('01_DECAY_TEST', 2):84.88250417,
                       ('01_DECAY_TEST', 3):67.41507249306659,
                       ('01_DECAY_TEST', 4):30.99138865654723,
                       ('01_DECAY_TEST', 5):0.0,
                       ('01_DECAY_TEST', 6):0.0,
                       ('01_DECAY_TEST', 7):0.0,
                       }






    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_out_activity_pipe_uniform_case_8():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    #test_case.print_all_nodes()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('branches', 1):14.60315678,
                       ('branches', 2):173.1009057,
                       ('branches', 3):139.8454576,
                       ('branches', 4):165.1433933,
                       ('branches', 5):157.5516907,
                       ('branches', 6):91.27385895,
                       ('branches', 7):125.6964107,
                       ('branches', 8):524.5331793,
                       ('branches', 9):145.4944061,
                       ('branches',10):434.2701606,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_out_activity_pipe_uniform_case_8_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
            "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    #test_case.print_all_nodes()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')

    expected_values = {('branches', 1):14.60315678,
                       ('branches', 2):173.1009057,
                       ('branches', 3):139.8454576,
                       ('branches', 4):165.1433933,
                       ('branches', 5):157.5516907,
                       ('branches', 6):91.27385895,
                       ('branches', 7):125.6964107,
                       ('branches', 8):524.5331793,
                       ('branches', 9):145.4944061,
                       ('branches',10):434.2701606,
                       }


    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_average_activity_pipe_uniform_case_8():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')

    expected_values = {('branches', 1):12.30771601,
                       ('branches', 2):94.06336437,
                       ('branches', 3):155.8824122,
                       ('branches', 4):169.0909436,
                       ('branches', 5):161.3177707,
                       ('branches', 6):127.8527093,
                       ('branches', 7):148.1366703,
                       ('branches', 8):561.4214839,
                       ('branches', 9):148.1899336,
                       ('branches',10):446.0595416,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

def test_average_activity_pipe_uniform_case_8_det():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
             "numerical_method" : "deterministic",
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')

    expected_values = {('branches', 1):12.30771601,
                       ('branches', 2):94.06336437,
                       ('branches', 3):155.8824122,
                       ('branches', 4):169.0909436,
                       ('branches', 5):161.3177707,
                       ('branches', 6):127.8527093,
                       ('branches', 7):148.1366703,
                       ('branches', 8):561.4214839,
                       ('branches', 9):148.1899336,
                       ('branches',10):446.0595416,
                       }

    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        print (value)
        print (test_dict[key])
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)

    assert all(check_vec)

#def test_out_activity_pipe_uniform_case_8_100s():
#    """
#    This function check the outlet vol activity in each pipe against calculated
#    ones
#    """
#    precision = 2.5e-2
#    test_setting = {
#            "steady_state": False,
#            "pipe_time_default": "uniform",
#            "mc_error_max": precision,
#            "t_bins":"100",
#                   }
#    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
#    test_case = FlunedSlCase(test_setting, test_path)
#    test_case.mc_solve()
#    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
#
#    expected_values = {('branches', 1):14.60315678,
#                       ('branches', 2):173.1009057,
#                       ('branches', 3):139.8454576,
#                       ('branches', 4):165.1433933,
#                       ('branches', 5):157.5516907,
#                       ('branches', 6):91.27385895,
#                       ('branches', 7):125.6964107,
#                       ('branches', 8):524.5331793,
#                       ('branches', 9):145.4944061,
#                       ('branches',10):434.2701606,
#                       }
#
#
#    a_vec = []
#    b_vec = []
#    check_vec = []
#
#    for key,value in expected_values.items():
#        print (value)
#        print (test_dict[key])
#        a_vec.append(value)
#        b_vec.append(test_dict[key])
#        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
#        check_vec.append(check)
#        if not check:
#            print('key: ', key)
#            print('expected value: ', value)
#            print('calculated value: ', test_dict[key])
#            print('check: ', check)
#
#    assert all(check_vec)
#
#def test_average_activity_pipe_uniform_case_8_100s():
#    """
#    This function check the outlet vol activity in each pipe against calculated
#    ones
#    """
#    precision = 5e-2
#    test_setting = {
#            "steady_state": False,
#            "pipe_time_default": "uniform",
#            "mc_error_max": precision,
#            "t_bins":"100",
#                   }
#    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
#    test_case = FlunedSlCase(test_setting, test_path)
#    test_case.mc_solve()
#    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
#
#    expected_values = {('branches', 1):12.30771601,
#                       ('branches', 2):94.06336437,
#                       ('branches', 3):155.8824122,
#                       ('branches', 4):169.0909436,
#                       ('branches', 5):161.3177707,
#                       ('branches', 6):127.8527093,
#                       ('branches', 7):148.1366703,
#                       ('branches', 8):561.4214839,
#                       ('branches', 9):148.1899336,
#                       ('branches',10):446.0595416,
#                       }
#
#    a_vec = []
#    b_vec = []
#    check_vec = []
#
#    for key,value in expected_values.items():
#        print (value)
#        print (test_dict[key])
#        a_vec.append(value)
#        b_vec.append(test_dict[key])
#        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
#        check_vec.append(check)
#        if not check:
#            print('key: ', key)
#            print('expected value: ', value)
#            print('calculated value: ', test_dict[key])
#            print('check: ', check)
#
#    assert all(check_vec)
#
#def test_out_activity_pipe_uniform_case_8_01s():
#    """
#    This function check the outlet vol activity in each pipe against calculated
#    ones
#    """
#    precision = 1e-2
#    test_setting = {
#            "steady_state": False,
#            "pipe_time_default": "uniform",
#            "mc_error_max": precision,
#            "t_bins":"0.1",
#                   }
#    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
#    test_case = FlunedSlCase(test_setting, test_path)
#    test_case.mc_solve()
#    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
#
#    expected_values = {('branches', 1):2.902337221,
#                       ('branches', 2):96.74457403,
#                       ('branches', 3):0.0,
#                       ('branches', 4):0.0,
#                       ('branches', 5):0.0,
#                       ('branches', 6):0.0,
#                       ('branches', 7):0.0,
#                       ('branches', 8):0.0,
#                       ('branches', 9):0.0,
#                       ('branches',10):0.0,
#                       }
#
#
#    a_vec = []
#    b_vec = []
#    check_vec = []
#
#    for key,value in expected_values.items():
#        print (value)
#        print (test_dict[key])
#        a_vec.append(value)
#        b_vec.append(test_dict[key])
#        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
#        check_vec.append(check)
#        if not check:
#            print('key: ', key)
#            print('expected value: ', value)
#            print('calculated value: ', test_dict[key])
#            print('check: ', check)
#
#    assert all(check_vec)
#
#def test_average_activity_pipe_uniform_case_8_01s():
#    """
#    This function check the outlet vol activity in each pipe against calculated
#    ones
#    """
#    precision = 1e-2
#    test_setting = {
#            "steady_state": False,
#            "pipe_time_default": "uniform",
#            "mc_error_max": precision,
#            "t_bins":"0.1",
#                   }
#    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_8)
#    test_case = FlunedSlCase(test_setting, test_path)
#    test_case.mc_solve()
#    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
#    # the values of all nodes except 1 and 8 have been simulated as it is not
#    # possible to calculate them analytically. The important is that the
#    # volumetric reaction rate functions correctly by passing a small fraction
#    # of the activity to the following node for every time step
#
#    expected_values = {('branches', 1):6.92957003800000,
#                       ('branches', 2):30.318102369740448,
#                       ('branches', 3):2.2020203975212866,
#                       ('branches', 4):9.969301804626852,
#                       ('branches', 5):0.0,
#                       ('branches', 6):0.725355094439575,
#                       ('branches', 7):1.4727141166538467,
#                       ('branches', 8):43.182804170000,
#                       ('branches', 9):0.00000000000000,
#                       ('branches',10):0.00000000000000,
#                       }
#
#
#
#    a_vec = []
#    b_vec = []
#    check_vec = []
#
#    for key,value in expected_values.items():
#        print (value)
#        print (test_dict[key])
#        a_vec.append(value)
#        b_vec.append(test_dict[key])
#        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
#        check_vec.append(check)
#        if not check:
#            print('key: ', key)
#            print('expected value: ', value)
#            print('calculated value: ', test_dict[key])
#            print('check: ', check)
#
#    assert all(check_vec)
#
def test_out_activity_pipe_uniform_case_9():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_9)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):2.990946033,
                       ('01_DECAY_TEST', 2):2.275720916,
                       ('01_DECAY_TEST', 3):169.7907339,
                       ('01_DECAY_TEST', 4):79.47320355,
                       ('01_DECAY_TEST', 5):8.149686861,
                       ('01_DECAY_TEST', 6):3.814588157,
                       ('01_DECAY_TEST', 7):1.785477535,
                       }



    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

def test_average_activity_pipe_uniform_case_9():
    """
    This function check the outlet vol activity in each pipe against calculated
    ones
    """
    precision = 1e-2
    test_setting = {
            "pipe_time_default": "uniform",
            "mc_error_max": precision,
                   }
    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_9)
    test_case = FlunedSlCase(test_setting, test_path)
    test_case.mc_solve()
    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
    expected_values = {('01_DECAY_TEST', 1):1.503041457,
                       ('01_DECAY_TEST', 2):2.617064889,
                       ('01_DECAY_TEST', 3):88.72024805,
                       ('01_DECAY_TEST', 4):118.9724181,
                       ('01_DECAY_TEST', 5):31.31740209,
                       ('01_DECAY_TEST', 6):5.710488022,
                       ('01_DECAY_TEST', 7):2.672883063,
                       }





    a_vec = []
    b_vec = []
    check_vec = []

    for key,value in expected_values.items():
        a_vec.append(value)
        b_vec.append(test_dict[key])
        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
        check_vec.append(check)
        if not check:
            print('key: ', key)
            print('expected value: ', value)
            print('calculated value: ', test_dict[key])
            print('check: ', check)


    assert all(check_vec)

#def test_out_activity_pipe_uniform_case_9_100s():
#    """
#    This function check the outlet vol activity in each pipe against calculated
#    ones
#    """
#    precision = 1e-2
#    test_setting = {
#            "steady_state": False,
#            "pipe_time_default": "uniform",
#            "mc_error_max": precision,
#            "t_bins":"100",
#                    }
#    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_9)
#    test_case = FlunedSlCase(test_setting, test_path)
#    test_case.mc_solve()
#    test_dict = test_case.get_all_attr('mc_out_activity_bq_m3')
#    expected_values = {('01_DECAY_TEST', 1):2.990946033,
#                       ('01_DECAY_TEST', 2):2.275720916,
#                       ('01_DECAY_TEST', 3):169.7907339,
#                       ('01_DECAY_TEST', 4):79.47320355,
#                       ('01_DECAY_TEST', 5):8.149686861,
#                       ('01_DECAY_TEST', 6):3.814588157,
#                       ('01_DECAY_TEST', 7):1.785477535,
#                       }
#
#
#
#    a_vec = []
#    b_vec = []
#    check_vec = []
#
#    for key,value in expected_values.items():
#        a_vec.append(value)
#        b_vec.append(test_dict[key])
#        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
#        check_vec.append(check)
#        if not check:
#            print('key: ', key)
#            print('expected value: ', value)
#            print('calculated value: ', test_dict[key])
#            print('check: ', check)
#
#
#    assert all(check_vec)
#
#def test_average_activity_pipe_uniform_case_9_100s():
#    """
#    This function check the outlet vol activity in each pipe against calculated
#    ones
#    """
#    precision = 1e-2
#    test_setting = {
#            "steady_state": False,
#            "pipe_time_default": "uniform",
#            "mc_error_max": precision,
#            "t_bins":"100",
#                    }
#    test_path = os.path.join(os.getcwd(),TEST_REL_PATH_CASE_9)
#    test_case = FlunedSlCase(test_setting, test_path)
#    test_case.mc_solve()
#    test_dict = test_case.get_all_attr('mc_average_activity_bq_m3')
#    expected_values = {('01_DECAY_TEST', 1):1.503041457,
#                       ('01_DECAY_TEST', 2):2.617064889,
#                       ('01_DECAY_TEST', 3):88.72024805,
#                       ('01_DECAY_TEST', 4):118.9724181,
#                       ('01_DECAY_TEST', 5):31.31740209,
#                       ('01_DECAY_TEST', 6):5.710488022,
#                       ('01_DECAY_TEST', 7):2.672883063,
#                       }
#
#
#
#
#
#    a_vec = []
#    b_vec = []
#    check_vec = []
#
#    for key,value in expected_values.items():
#        print (value)
#        print (test_dict[key])
#        a_vec.append(value)
#        b_vec.append(test_dict[key])
#        check =   isclose(value, test_dict[key], rtol=5*precision, atol=1e-08)
#        check_vec.append(check)
#        if not check:
#            print('key: ', key)
#            print('expected value: ', value)
#            print('calculated value: ', test_dict[key])
#            print('check: ', check)
#
#
#    assert all(check_vec)

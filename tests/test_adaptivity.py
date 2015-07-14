from test_beam_bend import create_problem, check_soln
from elastic import adaptive_execute


def test_beam_bend():
    import logging
    logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

    es, params = create_problem(0)
    params['error_threshold'] = 0.0000001
    params['max_iters'] = 15
    params['refine_fraction'] = 0.3
    params['dense'] = True
    result = adaptive_execute(2, es, params)
    check_soln(result)

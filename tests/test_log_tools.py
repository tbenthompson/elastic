from elastic.log_tools import log_exceptions, log_elapsed_time
from StringIO import StringIO
import logging

def capture_logs():
    str_logger = logging.getLogger()
    str_logger.setLevel(0)
    str_buffer = StringIO()
    handler = logging.StreamHandler(str_buffer)
    str_logger.addHandler(handler)
    return str_logger, str_buffer

def extract_log_text(str_buffer):
    str_buffer.flush()
    log_text = str_buffer.getvalue()
    return log_text

def test_log_exceptions():
    str_logger, str_buffer = capture_logs()

    @log_exceptions(str_logger)
    def function():
        raise Exception('WHOA')

    try:
        function()
    except:
        pass

    assert('Exception: WHOA' in extract_log_text(str_buffer))

def timing_test(do_timing):
    str_logger, str_buffer = capture_logs()

    @log_elapsed_time(str_logger, 'HI')
    def function(self):
        return 5

    class FakeTimer(object):
        def __init__(self):
            self.t = 0
        def time(self):
            out = self.t
            self.t += 1
            return out
    class A: pass

    a = A()
    a.params = dict()
    a.params['timing'] = do_timing
    a.params['timer'] = FakeTimer()
    result = function(a)
    assert(result == 5)

    log_text = extract_log_text(str_buffer)
    return log_text

def test_dont_log_elapsed_time():
    log_text = timing_test(False)
    assert(log_text == 'Starting HI\nFinished HI\n')

def test_log_elapsed_time():
    log_text = timing_test(True)
    assert(log_text == 'Starting HI\nFinished HI. Took 1\n')

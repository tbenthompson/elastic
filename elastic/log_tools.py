import traceback

def log_exceptions(logger):
    def wrapper(f):
        def wrapped(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                upper = ''.join(traceback.format_list(traceback.extract_stack()[:-1]))
                lower = '\n'.join(traceback.format_exc().split('\n')[1:-1])
                logger.error('Exception:\n' + upper + lower)
                raise
        return wrapped
    return wrapper

def log_elapsed_time(logger, text):
    if type(text) is str:
        text_fnc = lambda self, args: text
    else:
        text_fnc = text

    def wrapper(f):
        def wrapped(self, *args):
            base_msg = text_fnc(self, args)
            start_msg = 'Starting ' + base_msg
            logger.info(start_msg)
            if self.params['timing']:
                start_time = self.params['timer'].time()
            result = f(self, *args)
            finish_msg = 'Finished ' + base_msg
            if self.params['timing']:
                end_time = self.params['timer'].time()
                elapsed_time = end_time - start_time
                finish_msg += '. Took ' + str(elapsed_time)
            logger.info(finish_msg)
            return result
        return wrapped
    return wrapper


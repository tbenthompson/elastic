import traceback
import time

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

def log_elapsed_time(logger, text, timer = time):
    if type(text) is str:
        def text_fnc(*args, **kwargs):
            return text
    else:
        text_fnc = text

    def wrapper(f):
        def wrapped(*args, **kwargs):
            base_msg = text_fnc(*args, **kwargs)
            start_msg = 'Starting ' + base_msg
            logger.info(start_msg)
            start_time = timer.time()
            result = f(*args, **kwargs)
            finish_msg = 'Finished ' + base_msg
            end_time = timer.time()
            elapsed_time = end_time - start_time
            finish_msg += '. Took ' + str(elapsed_time)
            logger.info(finish_msg)
            return result
        return wrapped
    return wrapper


#TODO: Replace with the standard python logging utilities
class Logger(object):
    def output(self, data):
        pass

    def linear_solve_start(self, n_rows):
        self.output('Solving linear system with ' + str(n_rows) + ' rows.')

    def linear_solve_end(self):
        self.output('Finished solving linear system')

class StdoutLogger(Logger):
    def output(self, data):
        print(data)

#STUB
class FileLogger(Logger):
    def __init__(self, filename):
        pass

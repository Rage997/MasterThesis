import logging
logging.basicConfig(filename='example.log',
                format='%(levelname)s:%(message)s',
                filemode='w', level=logging.DEBUG, force=True)
logging.debug('This message should go to the log file')
logging.info('So should this')
logging.warning('And this, too')


print(logging.root.level)

print(logging.DEBUG)
#/usr/bin/env python3
"""
    Commonly used utilities

    Function    
    ---------------
    throw_error(msg):
        throw error and exit
    
    write_log(msg, lvl=20):
        write logging log to log file
    
    parse_tswildcard(tgt_time, wildcard):
        parse string with timestamp wildcard 
        to datetime object

"""
# ---imports---
import logging


# ---Module regime consts and variables---
print_prefix='lib.utils>>'

def parse_kwargs(kwargs, key, rtype='str'):
    '''
    parse kwargs from cfg string,
    e.g. 'a:1,3,5|b:2,c,d|c:3'
    '''
    if not(key in kwargs):
        return [] 
    
    kv_list=kwargs.split('|')
    
    for kv in kv_list: # a:1,3,5
        if key in kv:
            k,v=kv.split(':')
            vlist=v.split(',')
            if rtype=='str':
                return v
            elif rtype=='strlist':
                return vlist
            elif rtype=='int':
                return [int(itm) for itm in vlist]
            elif rtype=='float':
                return [float(itm) for itm in vlist]
            elif rtype=='bool':
                return [v.lower()=='true' for itm in vlist]
            
# ---Classes and Functions---
def throw_error(msg):
    '''
    throw error and exit
    '''
    logging.error(msg)
    exit()

def write_log(msg, lvl=20):
    '''
    write logging log to log file
    level code:
        CRITICAL    50
        ERROR   40
        WARNING 30
        INFO    20
        DEBUG   10
        NOTSET  0
    '''

    logging.log(lvl, msg)

def parse_tswildcard(tgt_time, wildcard):
    '''
    parse string with timestamp wildcard to datetime object
    '''
    seg_str=wildcard.split('@')
    parsed_str=''
    for seg in seg_str:
        if seg.startswith('%'):
            parsed_str+=tgt_time.strftime(seg)
        else:
            parsed_str+=seg
    return parsed_str


# ---Unit test---
if __name__ == '__main__':
    pass


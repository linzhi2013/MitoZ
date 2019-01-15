#!/usr/bin/env python3
import os
import sys
import re


def config2CL(config_file=None):
    para = ''
    start_flag = '[task={task}]'
    end_flag = '[task_{task}_input_end]'
    tmp = 0
    with open(config_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue
            if i.startswith('#'):
                continue

            if re.search(r'^task\s*=\s*(\w+)', i):
                m = re.search(r'^task\s*=\s*(\w+)')
                task = m.group(1)
                start_flag = start_flag.format(task=task)
                end_flag = end_flag.format(task=task)
                para = task
                continue

            if re.search(start_flag, i):
                tmp = 1
                continue

            if re.search(end_flag, i):
                tmp = 0
                break

            if tmp:
                m = re.search(r'^(\w+)\s*=\s*(.*)$', i)
                key = m.group(1)
                val = m.group(2)
                if val:
                    para += ' --{key} {val} '.format(key=key, val=val)

    return para

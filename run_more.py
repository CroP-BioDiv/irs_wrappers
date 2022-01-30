#!/usr/bin/env python3

import os.path
import importlib


def run(wrapper, filenames, output_results):
    if wrapper.endswith('.py'):
        wrapper = wrapper[:-3]
    module = importlib.import_module(wrapper)
    identify = module.identify

    results = []
    for filename in filenames:
        ret = identify(filename)
        print(filename, ret)
        ident = os.path.basename(filename).split('.', 1)[0]
        res = list(ret[0]) + list(ret[1]) if ret else [''] * 4
        results.append([ident] + res)

    if output_results:
        with open(output_results, 'w') as _out:
            _out.write('ident,ira_start,ira_end,irb_start,irb_end\n')
            for row in results:
                _out.write(f"{','.join(map(str, row))}\n")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run wrapper script with given sequences")
    parser.add_argument('wrapper', help='Wrapper script')
    parser.add_argument('-o', '--output-results', help='Output results in given CSV file')
    parser.add_argument('filenames', nargs='+', help='Sequence filenames')
    params = parser.parse_args()
    run(params.wrapper, params.filenames, params.output_results)

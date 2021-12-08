#!/usr/bin/python
import argparse
import sys
import glob
import multiprocessing as mp
from itertools import repeat
from pprint import pprint

import simplejson as json
from sqlalchemy import create_engine, Integer, Table, Column, MetaData, String
from sqlalchemy.dialects.postgresql import JSONB, insert
from sqlalchemy.exc import OperationalError

engine = create_engine(
    'postgresql+psycopg2://compound:asdf@lc-binbase-dev.czbqhgrlaqbf.us-west-2.rds.amazonaws.com/cts-test')
metadata = MetaData()
tables = {
    'pubchemc': Table('pubchemc', metadata, autoload_with=engine),
    'pubchems': Table('pubchems', metadata, autoload_with=engine)
}


def call_db(items, table, column, filename):
    pname = mp.current_process().name
    with engine.connect() as conn:
        insert_stmt = insert(tables[table]).values(items)
        do_nothing = insert_stmt.on_conflict_do_nothing(constraint=table + '_pk')

        try:
            with conn.begin():
                r = conn.execute(do_nothing)
        except OperationalError:
            try:
                with conn.begin():
                    r = conn.execute(do_nothing)
            except OperationalError as ex:
                with open('error.txt', 'a') as error:
                    error.write(filename + '\n')
                    error.flush()
                print(f'[{pname}]', ex)


def upload_chemspider():
    pass


def upload_pubchem(filename, folder):
    if filename.lower().startswith('compound'):
        id_label = 'PUBCHEM_COMPOUND_CID'
        table = {'name': 'pubchemc', 'column': 'compound'}
    else:
        id_label = 'PUBCHEM_SUBSTANCE_ID'
        table = {'name': 'pubchems', 'column': 'substance'}

    data = []
    with open(folder + '/' + filename, 'r') as fin:
        print(f'loading {folder}/{filename}')
        for f in fin:
            j = json.loads(f.strip())
            data.append({'id': j[id_label], 'file': filename, table['column']: j})

    print(f'grouping data in {filename}')
    dd = [data[i:i + 100] for i in range(0, len(data), 100)]
    # dd = [{'table': table['name'], 'column': table['column'], 'values': data[i:i + 10]} for i in
    #       range(0, len(data), 10)]

    p2 = mp.get_context('spawn').Pool(10)
    p2.map(call_db, zip(dd, table['name'], table['column'], repeat(f'{folder}/{filename}')), 1)


def main(params):
    if params['t']:
        last = 1
    else:
        last = -1

    files = [f.split('/')[-1] for f in glob.glob(params['input'] + '/*.json')[:last]]
    files.sort()

    with open(params['input']+'/done_cmp_json.txt', 'r') as done:
        skip = [s.strip() for s in done.readlines()]

    filtered = [f for f in files if f not in skip]

    for f in filtered:
        upload_pubchem(f, params['input'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Upload Pubchem .json files to RDS')
    parser.add_argument('-i', '--input', help='Folder with .json files', required=True)
    parser.add_argument('-t', action='store_true', help='Run in test mode.')

    try:
        args = parser.parse_args()
        args.input = args.input.rstrip('/')
        main(vars(args))
    except Exception as e:
        pprint(e.__cause__)
    finally:
        parser.print_help()

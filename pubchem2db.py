#!/usr/bin/python
import sys
import glob
import multiprocessing as mp
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


def call_db(item):
    pname = mp.current_process().name
    with engine.connect() as conn:
        insert_stmt = insert(tables[item['table']]).values(item['values'])
        do_nothing = insert_stmt.on_conflict_do_nothing(constraint=f"{item['table']}_pk")

        try:
            with conn.begin():
                r = conn.execute(do_nothing)
                # print(pname, r.inserted_primary_key_rows)
        except OperationalError:
            try:
                with conn.begin():
                    r = conn.execute(do_nothing)
                    # print(pname, r.inserted_primary_key_rows, '(retry)')
            except OperationalError as ex:
                with open('retry.sql', 'a') as retry:
                    retry.write(f'[{pname}] {do_nothing} ({do_nothing.mappings})\n')
                    retry.flush()
                print(f'[{pname}]', ex)


def upload(filename):
    # pname = mp.current_process().name
    fname = filename.split('/')[-1].strip()

    if fname.lower().startswith('compound'):
        id_label = 'PUBCHEM_COMPOUND_CID'
        table = {'name': 'pubchemc', 'column': 'compound'}
    else:
        id_label = 'PUBCHEM_SUBSTANCE_ID'
        table = {'name': 'pubchems', 'column': 'substance'}

    data = []
    with open(filename, 'r') as fin:
        print(f'loading {filename}')
        for f in fin:
            j = json.loads(f.strip())
            data.append({'id': j[id_label], 'file': fname, table['column']: f.strip()})

    print(f'grouping data in {filename}')
    dd = [{'table': table['name'], 'values': data[i:i + 10]} for i in
          range(0, len(data), 10)]

    p2 = mp.get_context('spawn').Pool(5)
    p2.map(call_db, dd, 1)


def main(args):
    last = -1
    if '-t' in args:
        last = 1

    files = glob.glob(f'{args[0]}*.json')
    files.sort()

    # with mp.get_context('spawn').Pool(1) as p:
    #     p.map(upload, files[:last], 1)
    for f in files[:last]:
        upload(f)


if __name__ == '__main__':
    main(sys.argv[1:])

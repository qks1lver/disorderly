# Disorderly.py
#### Compare protein sequences by their lengths and compositions
###### MIT License.
###### Requires Python 3+

To see the commands:
```
$ python3 disorderly.py -h
```

## How to use it?

### 1. Prepare your query
Put your query sequences in FASTA format and put them in a file

### 2. Prepare your database
Your database is made of sequences that you want to compare against.
This is also in FASTA format, but we need to convert it to a .disorderdb
database so it can be used to search against.
Generate a .disorderdb file from your database using the following
command:

```
$ python3 disorderly.py -v -fb path/to/your_database.fasta
```

__-v__ Verbose flag

__-fb__ Database FASTA file

This will generate __your_database.fasta.disorderdb__ in the same
folder as __your_database.fasta__

### 3. Search
Each of your queries is compared only to sequences of the same length in
the database.
Once a same-length sequence is found, the Euclidean distance between the
compositions
of your query and the database sequence is computed. The output contains
all the same-length sequences sorted by the Euclidean distance (low to
high).

__This search is distributed over all the available CPUs!__
```
$ python3 disorderly.py -v -i path/to/query.fasta -db path/to/your_database.fasta.disorderdb
```
__-i__ Your query sequences in FASTA

__-db__ The converted .disorderdb database

This will generate a __.csv__ with the same name as your query with a
bit of additional stuff (i.e. for __query.fasta__, the result will be
__query_search-20180816190934-ABCD.csv__). The __-v__ verbose flag will
tell you where your result is, which will be in the same directory as
your query)

### Alternatively, you can run everything all at once:
```
$ python3 disorderly.py -v -i query.fasta -fb your_database.fasta
```
The previous step-by-step instruction is meant to help you understand what is really
going on.

## Reading the result
##### Open the .csv file with a text editor or Excel
The format is (sequence IDs are the FASTA headers):

Queries|Hits|Distances
---|---|---
query-seq-1|database-seq-9|0.000
query-seq-1|database-seq-5|0.135
query-seq-1|database-seq-14|0.246
query-seq-2|database-seq-3|0.000
query-seq-2|database-seq-75|0.321


## How to get it? (Install)
#### No wheel currently :( , so just:
##### 1. Download the .zip
##### 2. Unpack it wherever you want
##### 3. Find disorderly.py under src/ and run as described above

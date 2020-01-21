# refget-py

I made two demos:

- [demo](demo.ipynb): shows how to create a redis or local database, load up some sequences and collections, and use refget to retrieve a fasta file
- [advanced](advanced.ipynb): shows how to do sequence-level comparisons between two checksums, revealing structure of the relationship (*e.g.* ordering differences, naming mismatches, sequence subsets).



```
docker run -it --network "host" mongo
```
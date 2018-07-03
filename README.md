# Hidden-Markov-Model
## compile
```
g++ source.cpp -o hmm
```

## Usage:
```
usage:
hmm --learn <datafile> <Nstates> <paramsfile>
hmm --display <paramsfile>
hmm --predict <paramsfile> <datafile> <outfile>
```
## datafile format
text file with columns separated by commas. Missing data is represented as 'nan'. Multiple sequences are separated by a space
### example
```

```

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
nan  ,  0.112439838825  ,  0.711295208763  ,  -0.0342975987291 
 -0.517634525609  ,  1.12031678197  ,  -0.460721543041  ,  nan 
 nan  ,  0.106853189055  ,  1.24320367838  ,  -0.716046567035 
 0.209240565226  ,  -0.0105126577448  ,  -0.506596276044  ,  -0.141423849337 
 nan  ,  -0.773251129163  ,  3.21599194358  ,  -0.424367460702 
 3.04618403765  ,  0.150161844676  ,  0.227244426573  ,  nan 
 
 -0.58414734158  ,  0.58393576546  ,  -0.434297118139  ,  0.827520474988 
 -1.10944733046  ,  -1.16979752471  ,  nan  ,  -0.569854526962 
 nan  ,  -1.34221277267  ,  0.0364329439943  ,  4.12070124973 
 2.17931198203  ,  -0.331823029065  ,  nan  ,  1.20427706898 
 -1.61527544312  ,  0.690075156579  ,  -1.59565001262  ,  0.540605470362 
 1.21372311715  ,  1.353835221  ,  nan  ,  -0.648687270656 
 2.88968892325  ,  -0.213838182977  ,  0.918798031667  ,  2.07283213804 
 
 1.17434485965  ,  nan  ,  0.796369086985  ,  -0.190846253049 
 nan  ,  0.00479335963347  ,  0.316787659577  ,  nan
```
data consists of three sequences. First one is of length 6, second is of length 7 and the third one is of length 2.
## options description
### --learn
learns MAP estimates of the model parameters.
### --display
displays model parameters in text
### --predict
computes the maximum likelihood estimates of the hidden states based on the model parameters supplied.

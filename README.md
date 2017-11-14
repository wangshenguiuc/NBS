## Network-based stratification of tumor mutations (NBS2.0)

### About
This directory contains the code and resources of the following paper:

<i>"RDifferentially mutated networks reveal complex mechanisms for cancer heterogeneity
". Under review. </i>

1. NBS2.0 is a network-based tumor stratification algorithm. It takes a somatic mutation profile as input. It then cluster tumors by integrating molecular networks with this somatic mutation profile.
2. The default network is a recently published PPI network (cite: inbioMap). This software also supported user network.
3. Our experimental results are based on TCGA data. A collection of pre-processed TCGA somatic mutation profiles can be obtained here (TODO: where should I put these data? dropbox?)
3. Please contact Sheng Wang swang141@illinois.edu if you have issues using this software.

### Overview of the Model
We introduce NBS2 algorithm to address the contamination of frequently mutated genes in network-based patient stratification. Our solution to alleviate the noise from frequently mutated genes is only using selective key gene interactions rather than all human gene interactions in network propagation. 

# Step 1. Construct a small and restrictive gene network
Briefly, we selected gene interactions according to their mutually exclusive patterns in cancer because these genes are most likely to be evolved in cancer development. Since frequently mutated genes are unlikely to exhibit mutually exclusive patterns with their neighbors, their interactions are discarded by our method, thus diminishing the influence of frequently mutated genes on other genes. 

# Step 2. Network propogation based on the constructed network
Somatic mutations for each patient are represented as a binary (1,0) states on genes, in which a ‘1’ indicates a gene for which mutation (a single-nucleotide base change or the insertion or deletion of bases) has occurred in the tumor relative to germ line. For each patient, we used network propagation to smooth the influence of each mutation over its neighborhood in the small and restrictive network.

# Step 3. Patient stratification based on the 'network-smoothed' profile
We then clustered the resulting matrix of ‘network-smoothed’ profile into a predefined number of subtypes (k=2,3,...,6) (Fig. 2b). Specifically, a compact representation of each patient was first calculated by applying singular value decomposition to ‘network-smoothed’ profile. K-means++ was then used to cluster patients into subtypes based on those compact vectors. 

# Step 4. Detecting differentially mutated subnetworks
We then detected signature subnetworks by comparing the smoothed mutation profiles of different subtypes and validated these subnetworks using various of experimental datasets. 

For further details, see Online Methods of our paper. 

### Sub-directories
  - this root directory contains impelmentation of the rationale model used for the beer review data. ``rationale.py`` implements the independent selection version and ``rationale_dependent.py`` implements the sequential selection version. See the paper for details.
  - [example_rationales](example_rationales) contains rationales generated for the beer review data. 
  - [ubuntu](ubuntu) contains alternative implementation for the AskUbuntu data.
  - [medical](medical) contains alternative implementation for medical report classification. 

<br>

### Data
  - **Pre-processed TCGA somatic mutation profile:** We provide pre-processed TCGA somatic mutation profile of 22 cancer types at [here](TODO: where should I put it). This should be sufficient for reproducing our results.   
   - **Sample network:** We use inBioMap as the sample network in this paper. The original network can be downloaded here. The constructed small and restrictive network of 22 cancer types can be accessed [here] (TODO: link)
  
**Important Note:** all data is for research-purpose only.

<br>

### Code Usage

To run the code, you need matlab (> 2015 I used) installed. Older matlab verision may not have an implementation of K-means++. Next:
  1. Clone the NBS2.0 repo
  2. Run `matlab main.m --help` to see all running options

Example run of beer review data:
```
THEANO_FLAGS='device=gpu,floatX=float32'        # use GPU and 32-bit float
python rationale.py                             # independent selection version
      --embedding /path/to/vectors              # path to somatic mutation profile (required)
      --train reviews.aspect0.train.txt.gz      # path to training set (required)
      --dev reviews.aspect0.heldout.txt.gz      # path to development set (required)        
      --load_rationale annotations.json         # path to rationale annotation for testing (required)
      --aspect 0                                # which aspect (-1 means all aspects)
      --dump outputs.json                       # dump selected rationales and predictions
      --sparsity 0.0003 --coherent 2.0          # regularizations
```

<br>

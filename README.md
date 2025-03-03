# pips-insilico
Code and example datasets for multimodal feature extraction and enzyme activity prediction.

## A. Set up repository
Clone github repository to local directory: 
```
git clone git@github.com:chicanagram/pips-insilico.git
```

To set up the required environment for running the code with conda, navigate into the pips-silico directory and install the listed modules from requirements.txt:
```
cd pips-insilico
conda create --name pips-insilico
pip install -r requirements.txt
source activate pips-insilico
```
## B. Dataset
The GOh1052 mutagenesis dataset is provided as an example at the location: 
```
data > ml_prediction > Input > GOh1052.csv
```
## C. Feature extraction software dependencies (non-python)
While the multimodal activity prediction model can be trained using the dataset provided (see Section B) using the open-source Autogluon library, full extraction of the dataset features requires some other supporting software.
This includes 1) YASARA Structure, 2) FoldX 3) Tango and 4) Waltz.
  
#### YASARA Structure
* Obtain the license for YASARA Structure from: http://www.yasara.org/products.htm
* Download and unzip the Yasara folder, following the installation instructions provided. The final folder should contain the appropriate Yasara executable. 
* Note the location of the Yasara executable, and update it in the file **yasara.py**. i.e. change the line 
```
yasaradir = '/Applications/YASARA.app/Contents/yasara/'
```
#### FoldX
* Obtain the license for FoldX from: http://foldxsuite.crg.eu/
* Download the FoldX executable. Update its full path in the file **yasara.py**. i.e. change the line
```
foldx_abspath = yasaradir + '/foldx_2025/foldx_20251231_mac'
```
#### Tango and Waltz
* Obtain the licenses for Tango and Waltz from: https://switchlab.org/software/
* Download the respective executables and place them in the subfolder **feature_extraction > aggregation**. 

## D. Perform Model Training & Testing
### Train
To perform training (e.g. on the GOh1052.csv dataset), run the following command from terminal:
```
python train.py
``` 
Results will be by default logged in the folder: examples > Output. 

The dataset filename or location be specified with the flags '-f' and '--data_folder', respectively. E.g. suppose the new base filename is 'GOh1030', and the new relative location is '../datasets/', run the command:
```
python train.py -f GOh1030mut --data_folder ../datasets/
```
To view all the arguments that can be specified, run: 
```
python train.py -h
```
These include 'num_bag_folds' (default=8), 'num_stack_levels' (default=1), 'save_model' etc. By default, the trained models will be saved at: data > ml_prediction > trained_models > <filename base>.

### Test
To perform testing (e.g. on the GOh1052mut.csv dataset), run the following command:
```
python test.py
``` 
Results will be logged in the folder: **data > ml_prediction > Output**. 

The dataset name and location can be set from commandline in a similar way as for the train.py script. 

## E. Perform train/test cross-validation on random splits
### Full dataset
To perform k-fold cross validation, run the command: 
```
python train_test_random.py
``` 
Here, the number of random splits is set as 4, but it can be modified using the 'n_splits' flag.

### Reduced dataset
Cross validation can also be done on the reduced datasets, whereby data is filtered according to one of: 1) normalized SIFT score (sift), 2) ShanMS score (shanMS), 3) distance to the active site (distance), of the position mutated. To perform data reduction, use the 'filt-by' flag.
The results of model cross validation will similarly be saved in the data > ml_prediction > Output folder. 
To filter using SIFT scores (default thresholds: 0.10 < x < 0.45)
```
python train_test_random.py --filt-by sift
``` 
To filter using ShanMS scores (default thresholds: 0.10 < x <= 0.45)
```
python train_test_random.py --filt-by shanMS
``` 
To filter using distance scores (default thresholds: 20 < x <= 50)
```
python train_test_random.py --filt-by distance
```
The filtering thresholds can also be adjusted using the flags 'filt-sift', 'filt-shanms', 'filt-dist', respectively. 
For example, to adjust the SIFT threshold to 0.20 < x <= 0.60, use the command:
```
python train_test_random.py --filt-by sift --filt-sift 0.20,0.60
```
whereby the upper and lower bound thresholds are specified, separated by a comma. 

## F. Obtain SIFT, Shannon Entropy & distance scores to perform retrospective data reduction
Code is provided in the **data_reduction** subdirectory to obtain the normalized SIFT, ShanMS and distance scores. To obtain SIFT and ShanMS scores, the relevant MSA needs to be added to the subfolder data > msa. Here, the file 'GOh1052_msa.fasta' is used, and the first sequence is the reference sequence. 
Obtaining of the SIFT output from the webserver at https://sift.bii.a-star.edu.sg/ has been automated. To obtain the distance scores, we use YASARA -- the free version should suffice for this particular calculation.
To obtain all the scores used for data reduction, run:
```
python get_data_reduction_scores.py
```

## G. Obtain binding, stability and aggregation features for multimodal ML prediction
Code is provided in the **feature_extraction** subdirectory to obtain the binding, stability and aggregation features using YASARA, FoldX, Waltz and Tango respectively.
Installation of the required software dependencies is outlined in Section C. Example can be found in the **data > feature_extraction** subfolder showing outputs from running the YASARA and FoldX feature extraction pipelines on the first 5 mutants (sequentially) of GOh1052, as well as Tango and Waltz on the full set of mutants. 
```
python get_binding_stability_agg_features.py
```

## H. Feature importances for autogluon models
Finally, feature importances for the multimodal models may be obtained using Autogluon's built in module by running: 
```
python get_feature_importance.py
```
The model analysed may be specified by modify8ing the **'model_name'** variable the script. 
It is currently set to output feature importances for the 'WeightedEnsemble_L3' model saved in the **data > trained_models > GOh1052** subdirectory, the model trained with the best performance. 


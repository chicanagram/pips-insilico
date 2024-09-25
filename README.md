# pips-insilico
Code and example datasets for running the in-silico models developed during the PIPS project

## Set up repository
Clone github repository to local directory: 
```
git clone git@github.com:chicanagram/pips-insilico.git
```

To set up the required environment for running the code with conda, navigate into the pips-silico directory and install the listed modules from requirements.txt:
```
cd pips-insilico
conda create --name pips-insilico
pip install -r requirements.txt
```

## Example datasets
The GOh1052 mutagenesis dataset is provided as an example at the location: 
```
examples > Input > GOh1052mut.csv
```
This dataset has also been randomly split into a a train and test set, found at the same location
```
examples > Input > GOh1052mut_train.csv
examples > Input > GOh1052mut_test.csv
```
## Perform Model Training
To perform training (e.g. on the GOh1052mut_train.csv dataset), run the following command from terminal:
```
python train.py
``` 
Results will be by default logged in the folder: examples > Output. 

If there is a need to change the dataset filename or location, this can be specified with the flags '-f' and '--data_folder', respectively. E.g. suppose the new base filename is 'GOh1030mut', and the new relative location is '../datasets/', run the command:
```
python train.py -f GOh1030mut --data_folder ../datasets/
```
To view all the arguments that can be specified, run: 
```
python train.py -h
```
These include 'num_bag_folds' (default=8), 'num_stack_levels' (default=1), 'save_model' etc. By default, the trained models will be saved at: examples > trained_models > <filename base>.

##  Perform Model Testing
To perform testing (e.g. on the GOh1052mut_test.csv dataset), run the following command:
```
python train.py
``` 
Results will be logged in the folder: examples > Output. 

The dataset name and location can be set from commandline in a similar way as for the train.py script. 

## Perform k-fold train/test cross-validation
To perform k-fold cross validation, run the command: 
```
python train_test_kfold.py --n_splits 4
``` 
Here, the number of random splits is set as 4, but it can be modified. 
The arguments available for running the train.py script are apply here as well. 

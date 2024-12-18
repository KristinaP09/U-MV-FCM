# U-MV-FCM
Matlab code for the paper "Unsupervised Multiview Fuzzy C-Means Clustering Algorithm"

codes written by: Kristina P. Sinaga

### Method Description

Our method can find the optimal number of clusters with high-quality cluster centroids in an efficient manner, called U-MV-FCM. The proposed U-MV-FCM can automatically produce an optimal number of clusters and simultaneously improve the accuracy rate without a parameter setting.The view-points on  the multiview data are assigned as the initial cluster centers and the proposed U-MV-FCM can directly reduce the number of clusters and automatically produce an optimal number of clusters.

### Environment

U-MV-FCM was developed in MATLAB 2020a

### Usage

We provided demo on Syn500 for users. To run this demo, please load the script 'demo_URMVFCM.m' into your MATLAB programming environment and click 'run'.


### Parameters

There are three parameters in our method, i.e., 'eta_1', 'eta_2', "eta_3" and the exponent parameter 'beta'. These four parameters are estimated. However, users can change their value in 'demo_URMVFCM.m' .

### Input and output

Users can change the input file directory and output file directory by changing the input data matrix  in 'demo_URMVFCM.m', respectively.

### Citation

```

@article{hussain2023unsupervised,
  title={Unsupervised multiview fuzzy c-means clustering algorithm},
  author={Hussain, Ishtiaq and Sinaga, Kristina P and Yang, Miin-Shen},
  journal={Electronics},
  volume={12},
  number={21},
  pages={4467},
  year={2023},
  publisher={MDPI}
}

```

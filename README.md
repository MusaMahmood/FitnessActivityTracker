# Fitness Tracking Android Application Project
The previous classification system was experimental and did not work most of the time. Most of the code from the previous classifier was thrown out, and a new set of training data recorded. The classes were narrowed down to 3:
- 0 Idle, little to no movement
- 1 Idle, threshold exceeded
- 2 Walking 
- 3 Jogging or Running. 
With the aid of a threshold based system, we were able to attain a higher accuracy with an updated feature set. The total number of features was increased to 13 per input data set. As there are 8 total data sets (Acceleration on the X Y Z axes, and the resultant vector, as well as Gyroscope rotational data on the X Y Z axes and resultant vector), this gives us a total of 104 features. 
Training data was recorded in three sessions, one for general idle use, one recording of the user walking and another of the user jogging and sprinting. The total length of the recordings is about 21 minutes. The data is divided into 50-sample (or about 3 seconds at 16.6Hz) intervals and the features are extracted and recorded in a matrix, along with a corresponding class label. 

With the new dataset and updated features, we attained the following classifier accuracy prediction using MATLAB's Classification Learner:
![MATLAB Classification Learner Fine-kNN](https://cloud.githubusercontent.com/assets/25552144/25310240/057995c8-27ae-11e7-818f-248a2494c1d1.png)

The classification learner suggests a 99.6% accuracy, but in practice it is slightly lower than this. Using the 50% holdout validation option, we see the failure rate does not change, which is a good sign that our classification system is good.
![image](https://cloud.githubusercontent.com/assets/25552144/25310285/37f606c0-27af-11e7-857f-b0538b44a435.png)

Human testing shows that the accuracy of the classifier is good, but not perfect, and can definitely still be improved. 
You can see how the thresholds separate important features in the following plot of the accelerometer resultant:
![MATLAB Accelerometer Training Data](https://cloud.githubusercontent.com/assets/25552144/25310262/b1caca86-27ae-11e7-8b7d-f964c46bf8bb.png)

This can also be seen in the gyroscope data. 

![MATLAB Gyroscope Training Data](https://cloud.githubusercontent.com/assets/25552144/25310270/ee5ed046-27ae-11e7-9a60-9144b10401e3.png)

Currently, the app will record acceleration and rotation data, and will store it in a data buffer, and call the classificaiton function every 3 seconds. The classified data result is stored in a csv file for easy access. These data points are not currently timestamped, although that would make it far more useful for tracking activity over time periods shorter than 1 day (hour by hour progress). 

Test Plan:
The plan to test this app was to simply run the app in the background while keeping track of activity on paper, so that the tester could check the accuracy of the recorded classes against the expected classes. These simple tests were conducted, and provided a greater than 90% accuracy on all trials, which was the original aim. However, this was only tested on a single device (Samsung Galaxy S7), and does depend heavily on the accuracy of the integrated motion sensor chip.
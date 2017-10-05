package musamahmood.fitnessapp;

import android.content.Context;
import android.content.Intent;
import android.content.pm.ActivityInfo;
import android.content.pm.PackageManager;
import android.graphics.Color;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.media.MediaPlayer;
import android.net.Uri;
import android.os.Bundle;
import android.os.Environment;
import android.support.v4.app.ActivityCompat;
import android.support.v4.content.ContextCompat;
import android.support.v4.content.FileProvider;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.Display;
import android.view.View;
import android.view.WindowManager;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.ImageView;
import android.widget.TextView;
import android.widget.Toast;
import android.widget.ToggleButton;

import com.androidplot.Plot;
import com.androidplot.util.Redrawer;
import com.androidplot.xy.BoundaryMode;
import com.androidplot.xy.LineAndPointFormatter;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.XYPlot;
import com.androidplot.xy.XYStepMode;
import com.opencsv.CSVWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Locale;

// Android Plot
// Google Activity Recognition
// Google MediaPlayer Service

public class MainActivity extends AppCompatActivity implements SensorEventListener {

    //Initialize native library for C++ integration:
    static {
        System.loadLibrary("native-lib");
        System.loadLibrary("resultant-lib");
    }

    private static final int MY_PERMISSIONS_WRITE_STORAGE = 860;

    //Debug Tag
    private final static String TAG = MainActivity.class.getSimpleName();

    private TextView ActivityDetected_textview;
    private TextView mTextMessage;
    private ImageView image;
    public static final String PACKAGE_NAME = "vcucmsc355.fitnessapp;";

    // Sensor Variables for Raw Data Acquisition
    private SensorManager mSensorManager = null;
    private Sensor mAccelerometerSensor = null;
    private Sensor mGyroscopeSensor = null;
    private double[] mMotionSensorDataBuffer = new double[6];

    private double[] mAccXData = new double[50];
    private double[] mAccYData = new double[50];
    private double[] mAccZData = new double[50];
    private double[] mGyrXData = new double[50];
    private double[] mGyrYData = new double[50];
    private double[] mGyrZData = new double[50];

    private int countDataPoints = 0;

    private boolean mAccDataPresent = false;
    private boolean mGyroDataPresent = false;

    //Android Plot Series
    private SimpleXYSeries mAccResSeries = null; //Data for plotting accelerometer resultant
    private SimpleXYSeries mGyroResSeries = null; //Data for plotting gyroscope resultant
    private XYPlot mMotionPlot = null;
    public static final int DATA_HISTORY = 1000;

    //Redrawer for Android Plot plugin
    private Redrawer mRedrawer;

    //Data Acquisition Variables:
    private boolean mDataAcquisitionEnabled = false;

    //UI Variables
    private Button mActivityButton;
    private Button mStatisticsButton;
    private TextView mStepCounterText;
    private TextView mActivityTrackingStatus;

    // MediaPlayer Variable
    MediaPlayer happiness;
    MediaPlayer ringtone;
	MediaPlayer musicList;
	ArrayList<Integer> playlist;

    // Variables for Data Storage as CSV:
    private boolean recordRawData = false;
    private boolean fileSaveInitialized = false;
    private CSVWriter csvWriter;
    private CSVWriter csvWriter2;
    private File root;
    private File file;
    private String[] csvWriteData = new String[6];
    private String[] csvWriteData2 = new String[5];
    private String fileTimeStamp = "";

    //Variables for getting Data Rate
    private long mLastTime;
    private long mLastTime2;
    private long mCurrentTime;
    private long mCurrentTime2;
    private int points = 0;
    private double dataRate;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        setRequestedOrientation(ActivityInfo.SCREEN_ORIENTATION_PORTRAIT);
        // MediaPlayer
        happiness = MediaPlayer.create(this, R.raw.happiness);
        happiness.setLooping(true);
        ringtone = MediaPlayer.create(this, R.raw.ringtone);
        ringtone.setLooping(true);
		playlist = new ArrayList<>();
		playlist.add(R.raw.happiness);
		playlist.add(R.raw.ringtone);
		musicList = MediaPlayer.create(this,playlist.get(0));
		
        image = (ImageView) this.findViewById(R.id.imageView1) ;
        ActivityDetected_textview = (TextView) findViewById(R.id.ActivityDetected_textview);
            // Bottom Navigation
        int permissionCheck = ContextCompat.checkSelfPermission(MainActivity.this, android.Manifest.permission.WRITE_EXTERNAL_STORAGE);
        if(permissionCheck!= PackageManager.PERMISSION_GRANTED) {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{android.Manifest.permission.WRITE_EXTERNAL_STORAGE}, MY_PERMISSIONS_WRITE_STORAGE);
        }

        mTextMessage = (TextView) findViewById(R.id.message);
        mSensorManager = (SensorManager) getApplicationContext().getSystemService(Context.SENSOR_SERVICE);
        //Check for presence of accelerometer:
        if (mAccelerometerSensor==null) {
            if(mSensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER)!=null) {
                mAccelerometerSensor = mSensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER);
                Log.d(TAG, "Accelerometer Enabled");
            } else {
                Toast.makeText(this,"Accelerometer Not Found!",Toast.LENGTH_LONG).show();
                Log.d(TAG, "Accelerometer Not Present");
            }
        }
        countDataPoints = 0;
        //Check for gyroscope
        if(mGyroscopeSensor==null) {
            if(mSensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE)!=null) {
                mGyroscopeSensor = mSensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE);
                Log.d(TAG, "Gyro Sensor Enabled");
            } else {
                Toast.makeText(this,"Gyroscope Not Found!",Toast.LENGTH_LONG).show();
                Log.d(TAG, "Gyro Sensor Not Present");
            }
        }
        //Initialize Buttons:
        mActivityButton = (Button) findViewById(R.id.trackActivityButton);
        mStatisticsButton = (Button) findViewById(R.id.statisticsButton);
        mStepCounterText = (TextView) findViewById(R.id.stepCounterText);
        mActivityTrackingStatus = (TextView) findViewById(R.id.ActivityTrackingStatus);
        // Begin listening for sensors, delay set to "game" (fast) sampling rate.
        enableSensorTracking();
        mMotionPlot = (XYPlot) findViewById(R.id.motionPlot);
        mMotionPlot.setRangeBoundaries(0, 6, BoundaryMode.FIXED);
        mMotionPlot.setDomainBoundaries(0,DATA_HISTORY,BoundaryMode.FIXED);
        mMotionPlot.setRangeValueFormat(new DecimalFormat("#.#"));
        mMotionPlot.setDomainValueFormat(new DecimalFormat("#"));
        mMotionPlot.setRangeLabel("Acceleration in Gs and Gyro in rads");
        mAccResSeries = new SimpleXYSeries("Acc Resultant (gs)");
        mGyroResSeries = new SimpleXYSeries("Gyro Resultant (rad)");
        mAccResSeries.useImplicitXVals();
        mGyroResSeries.useImplicitXVals();
        LineAndPointFormatter accLine = new LineAndPointFormatter(Color.rgb(100, 100, 200), null, null, null);
        LineAndPointFormatter gyroLine = new LineAndPointFormatter(Color.parseColor("#000000"),null, null, null);
        mMotionPlot.addSeries(mAccResSeries, accLine);
        mMotionPlot.addSeries(mGyroResSeries, gyroLine);
        mMotionPlot.setDomainStepMode(XYStepMode.INCREMENT_BY_VAL);
        mMotionPlot.setDomainStepValue(DATA_HISTORY/5);
        mMotionPlot.setDomainLabel("Sample Index");
        //Get Screen Refresh Rate for Redrawer
        Display display = ((WindowManager)getSystemService(Context.WINDOW_SERVICE)).getDefaultDisplay();
        final int refreshRate = (int)display.getRefreshRate();
        Log.d(TAG, "Max Refresh Rate: "+String.valueOf(refreshRate));
        mRedrawer = new Redrawer(Arrays.asList(new Plot[]{mMotionPlot}), refreshRate ,false);
        mMotionPlot.getDomainLabelWidget().getLabelPaint().setColor(Color.BLACK);
        mMotionPlot.getDomainLabelWidget().getLabelPaint().setTextSize(32);
        mMotionPlot.getRangeLabelWidget().getLabelPaint().setColor(Color.BLACK);
        mMotionPlot.getRangeLabelWidget().getLabelPaint().setTextSize(32);
        mMotionPlot.getGraphWidget().getDomainTickLabelPaint().setColor(Color.BLACK);
        mMotionPlot.getGraphWidget().getRangeTickLabelPaint().setColor(Color.BLACK);
        mMotionPlot.getGraphWidget().getDomainTickLabelPaint().setTextSize(36);
        mMotionPlot.getGraphWidget().getRangeTickLabelPaint().setTextSize(36);
        mMotionPlot.getGraphWidget().getDomainGridLinePaint().setColor(Color.WHITE);
        mMotionPlot.getGraphWidget().getRangeGridLinePaint().setColor(Color.WHITE);
        mMotionPlot.getLegendWidget().getTextPaint().setColor(Color.BLACK);
        mMotionPlot.getLegendWidget().getTextPaint().setTextSize(32);
        mMotionPlot.getTitleWidget().getLabelPaint().setTextSize(32);
        mMotionPlot.getTitleWidget().getLabelPaint().setColor(Color.BLACK);
        Button exportButton = (Button) findViewById(R.id.exportButton);
        exportButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
                    saveDataFile(true);
                } catch (IOException e) {
                    Log.e(TAG, "IOException in saveDataFile");
                    e.printStackTrace();
                }
                Context context = getApplicationContext();
                Uri uri = FileProvider.getUriForFile(context,context.getApplicationContext().getPackageName()+".provider", file);//Uri.fromFile(file);
                Intent exportData = new Intent(Intent.ACTION_SEND);
                exportData.putExtra(Intent.EXTRA_SUBJECT,"Motion Sensor Data Details");
                exportData.putExtra(Intent.EXTRA_STREAM, uri);
                exportData.setType("text/html");
                startActivity(exportData);
            }
        });
        ToggleButton trackRawSensorActivityToggleButton = (ToggleButton) findViewById(R.id.recordRawDataToggleButton);
        trackRawSensorActivityToggleButton.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean b) {
                recordRawData = b;
            }
        });
        mActivityButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(mDataAcquisitionEnabled) {
                    //Disable activity tracking.
                    disableSensorTracking();
                } else {
                    enableSensorTracking();
                }
            }
        });
        mStatisticsButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                Intent intent = new Intent(getApplicationContext(), DataControlActivity.class);
                startActivity(intent);
            }
        });
//        Button button2 = (Button) findViewById(R.id.button2);
//        button2.setOnClickListener(new View.OnClickListener() {
//            @Override
//            public void onClick(View v) {
//                long AA = System.currentTimeMillis();
//                disableSensorTracking();
//                while (true) {
//                    if(AA+5200<System.currentTimeMillis()) {
//                        enableSensorTracking();
//                        break;
//                    }
//                }
//            }
//        });
        mLastTime = System.currentTimeMillis();
        mLastTime2 = System.currentTimeMillis();
        stringFromJNI();
    }

    public static String getDateStamp() {
        return new SimpleDateFormat("yyyy.MM.dd", Locale.US).format(new Date());
    }

    /**
     *
     * @param terminate: if true, closes file writers, and stops writing to drive.
     * @throws IOException
     */
    public void saveDataFile(boolean terminate) throws IOException {
        if(terminate && fileSaveInitialized) {
            csvWriter.flush();
            csvWriter.close();
            csvWriter2.flush();
            csvWriter2.close();
        }
    }

    /**
     * Method for creating data file. It checks for specific file name to see if it exists.
     * The format is 'MotionSensorData_yyyy.MM.dd.csv'. If it already exists, it appends
     * to the end of the file.
     * @throws IOException
     */
    public void saveDataFile() throws IOException {
        root = Environment.getExternalStorageDirectory();
        fileTimeStamp = "MotionSensorData_"+getDateStamp();
        if(root.canWrite()) {
            File directory = new File(root.getAbsolutePath() + "/MotionData");
            directory.mkdirs();
            file = new File(directory,fileTimeStamp+".csv");
            File file2 = new File(directory,fileTimeStamp+"Classified.csv");
            if(file.exists() &&!file.isDirectory()) {
                Log.d(TAG, "File "+file.toString()+" already exists - appending data");
                FileWriter fileWriter = new FileWriter(file, true);
                csvWriter = new CSVWriter(fileWriter);
            } else {
                csvWriter = new CSVWriter(new FileWriter(file));
            }
            if(file2.exists() && !file2.isDirectory()) {
                FileWriter fileWriter = new FileWriter(file2, true);
                csvWriter2 = new CSVWriter(fileWriter);
            } else {
                csvWriter2 = new CSVWriter(new FileWriter(file2));
            }
        }
        fileSaveInitialized = (csvWriter!=null) && (csvWriter2!=null);
    }

    /**
     *
     * @param motionSensorData array to be saved to file.
     */
    public void saveDataFile(double[] motionSensorData) {
        if(fileSaveInitialized) {
            for (int i = 0; i < motionSensorData.length; i++) {
                csvWriteData[i] = motionSensorData[i]+"";
            }
            csvWriter.writeNext(csvWriteData, false);
        } else {
            try {
                saveDataFile();
            } catch (IOException e) {
                Log.e(TAG, "IOException in saveDataFile()");
                e.printStackTrace();
            }
        }
    }

    /**
     * @param y classified data array to be written to file.
     */
    public void saveClassifiedData(double y) {
        if(fileSaveInitialized) {
            String[] csvWriteData = new String[1];
            csvWriteData[0] = y+"";
            csvWriter2.writeNext(csvWriteData,false);
        } else {
            try {
                saveDataFile();
            } catch (IOException e) {
                Log.e(TAG, "IOException in saveClassifiedData()");
                e.printStackTrace();
            }
        }
    }

    @Override
    protected void onStart() {
        super.onStart();
//        mGoogleApiClient.connect();
    }

    @Override
    protected void onResume() {
        //Create File:
        if(!fileSaveInitialized) {
            try {
                saveDataFile();
            } catch (IOException e) {
                Log.e(TAG, "IOException in saveDataFile");
                e.printStackTrace();
            }
        }
        //Test JNI Function Call (Temporary):
        double Y = jniActivityTrackerInit();
        Log.d(TAG,"Initialize JNI Tracker Result: ["+String.valueOf(Y)+"]");
        super.onResume();
        if(!mDataAcquisitionEnabled)
            enableSensorTracking();
        mRedrawer.start();
    }

    @Override
    protected void onPause() {
//        LocalBroadcastManager.getInstance(this).unregisterReceiver(ADBroadcastReceiver);
//        fileSaveInitialized=false;
        mRedrawer.pause();
        super.onPause();
    }

    @Override
    protected void onDestroy() {
        try {
            saveDataFile(true);
        } catch (IOException e) {
            Log.e(TAG, "IOException in saveDataFile(terminate)");
            e.printStackTrace();
        }
        if(mDataAcquisitionEnabled)
            disableSensorTracking();
        mRedrawer.finish();
        super.onDestroy();
    }

    /**
     * Internal callback from SensorManager when new sample of an enabled sensor is detected.
     * @param sensorEvent contains sensor data: sensor ID, timestamp, accuracy, and measured values
     */
    @Override
    public synchronized void onSensorChanged(SensorEvent sensorEvent) {
        if(sensorEvent.sensor==mAccelerometerSensor) {
            mAccDataPresent = true;
            double accArray[] = new double[sensorEvent.values.length];
            for (int i = 0; i < sensorEvent.values.length; i++) {
                accArray[i] = sensorEvent.values[i];
                mMotionSensorDataBuffer[i] = sensorEvent.values[i];
            }
            //Converts from m/s^2 to Gs
            double accRes = jniResultantAcc(accArray);
            if(mAccResSeries.size() > DATA_HISTORY) {
                mAccResSeries.removeFirst();
            }
            mAccResSeries.addLast(null,accRes);
            getDataRateAcc();
        } else if (sensorEvent.sensor == mGyroscopeSensor) {
            mGyroDataPresent = true;
            double gyroArray[] = new double[sensorEvent.values.length];
            for (int i = 0; i < sensorEvent.values.length; i++) {
                gyroArray[i] = sensorEvent.values[i];
                mMotionSensorDataBuffer[i+3] = sensorEvent.values[i];
            }
            double gyroRes = jniResultantGyro(gyroArray);
            if(mGyroResSeries.size() > DATA_HISTORY) {
                mGyroResSeries.removeFirst();
            }
            mGyroResSeries.addLast(null,gyroRes);
            getDataRate();
        }
        //Once we receive both Acc/Gyro Data
        if(mAccDataPresent && mGyroDataPresent) {
            mAccDataPresent = false;
            mGyroDataPresent = false;
            addToDataBuffers(mMotionSensorDataBuffer);
            if(recordRawData)
                saveDataFile(mMotionSensorDataBuffer);//save combined array to file.
        }
    }

    /**
     * Adds motion sensor samples to buffer.
     * @param dataBuffer array of single motion data samples
     */
    private void addToDataBuffers(double[] dataBuffer) {
        //shift back:
        System.arraycopy(mAccXData,1,mAccXData,0,49);
        System.arraycopy(mAccYData,1,mAccYData,0,49);
        System.arraycopy(mAccZData,1,mAccZData,0,49);
        System.arraycopy(mGyrXData,1,mGyrXData,0,49);
        System.arraycopy(mGyrYData,1,mGyrYData,0,49);
        System.arraycopy(mGyrZData,1,mGyrZData,0,49);
        //add to end
        mAccXData[49] = dataBuffer[0];
        mAccYData[49] = dataBuffer[1];
        mAccZData[49] = dataBuffer[2];
        mGyrXData[49] = dataBuffer[3];
        mGyrYData[49] = dataBuffer[4];
        mGyrZData[49] = dataBuffer[5];
        countDataPoints++;
        if(countDataPoints==50) {
            countDataPoints=0;
            classifyData();

        }
    }

    /**
     * This function calls the classifier using the stored data.
     */
    private void classifyData() {
        final double Y = jniActivityClassifier(mAccXData,mAccYData,mAccZData,mGyrXData,mGyrYData,mGyrZData);
        Log.e(TAG,"jniActivityClassifier :: Y = "+String.valueOf(Y));
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                ActivityDetected_textview.setText("Current Activity: "+String.valueOf(Y));
            }
        });
        saveClassifiedData(Y);
    }

    /**
     * This method ends sensor data acquisition.
     */
    public void disableSensorTracking() {
        mSensorManager.unregisterListener(this);
        mDataAcquisitionEnabled = false;
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                mActivityButton.setText(getString(R.string.enable_sensor));
                mActivityTrackingStatus.setText(R.string.tracking_disabled);
            }
        });
        try {
            saveDataFile(true);
        } catch (IOException e) {
            Log.e(TAG, "IOException in saveDataFile(terminate)");
            e.printStackTrace();
        }
    }

    /**
     * This method enables sensor tracking. Currently the Accelerometer and Gyroscope are set
     * to sample at 16.6Hz or once every 60ms. Upon sensor detection,
     */
    public void enableSensorTracking() {
        mSensorManager.registerListener(this, mAccelerometerSensor, SensorManager.SENSOR_DELAY_UI); //~16.6Hz
        mSensorManager.registerListener(this, mGyroscopeSensor, 60000); //16.6Hz
        mDataAcquisitionEnabled = true;
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                mActivityButton.setText(getString(R.string.disable_sensor));
                mActivityTrackingStatus.setText(getString(R.string.tracking_enabled));
            }
        });
    }

    /**
     * Simple method for acquiring sensor data rate over a period of roughly 5 seconds.
     */
    private void getDataRate() {
        mCurrentTime = System.currentTimeMillis();
        points++;
        if(mCurrentTime > (mLastTime+5000)) {
            dataRate = points/5;
            points = 0;
            mLastTime = mCurrentTime;
            Log.i(TAG, "Gyroscope Data Rate: "+String.valueOf(dataRate)+"samples/second");
        }
    }

    private int points2 = 0;
    private void getDataRateAcc() {
        mCurrentTime2 = System.currentTimeMillis();
        points2++;
        if(mCurrentTime2 > (mLastTime2+5000)) {
            final double dataRateA = points2/5;
            points2 = 0;
            mLastTime2 = mCurrentTime2;
            Log.i(TAG, "Accelerometer Data Rate: "+String.valueOf(dataRateA)+"samples/second");
        }
    }

    /**
     * onAccuracyChanged(Sensor sensor, int accuracy)
     * Called when the accuracy of the registered sensor has changed.
     */
    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {
        // Don't use this
    }

    /*
     *****************************************
     *  C++ JNI Function Initializers:
     */

    /**
     * JNI Test Function
     * @return Test String from C++ Code
     */
    public native String stringFromJNI();

    /**
     * Simple functions for calculating resultant using JNI framework
     * @param array of length [3]
     * @return resultant of array data; Converts m/s^2 to Gs
     */
    public static native double jniResultantAcc(double[] array);

    /**
     * Simple functions for calculating resultant using JNI framework
     * @param array of length [3]
     * @return resultant of array data
     */
    public native double jniResultantGyro(double[] array);

    /**
     *
     * @param array Array of motion sensor data of length 1000
     * @param len length of the vector values being used. Here, the length of the array is static
     *            (1000), but only part of the array is relevant to the calculation.
     * @return mean of array from index '0' to index 'len'.
     */
    public static native double jniMean(double[] array, int len);

    /**
     * @param accX, accY, accZ - data from accelerometer
     * @param gyrX, gyrY ,gyrZ - data from gyroscope
     * @return Best fit of activity data using K-nearest neighbor algorithm
     */
    public native double jniActivityClassifier(double[] accX, double[] accY, double[] accZ,
                                              double[] gyrX, double[] gyrY, double[] gyrZ);

    /**
     *
     * @return vector 'Y' of length 5 for test call of activityTracker C++ method. It should return
     * all 1s.
     */
    public native double jniActivityTrackerInit();
}

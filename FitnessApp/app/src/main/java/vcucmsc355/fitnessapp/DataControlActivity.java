package vcucmsc355.fitnessapp;

import android.app.Activity;
import android.content.Context;
import android.graphics.Color;
import android.os.Bundle;
import android.os.Environment;
import android.util.Log;
import android.view.Display;
import android.view.WindowManager;
import android.widget.TabHost;

import com.androidplot.Plot;
import com.androidplot.util.Redrawer;
import com.androidplot.xy.BoundaryMode;
import com.androidplot.xy.LineAndPointFormatter;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.XYPlot;
import com.androidplot.xy.XYStepMode;
import com.opencsv.CSVReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

/**
 * Created by mahmoodms on 3/25/2017.
 */

public class  DataControlActivity extends Activity {
    private final static String TAG = DataControlActivity.class.getSimpleName();
    private final static int DATA_HISTORY = 1440; //Mins/24Hrs
    private CSVReader csvReaderRawData;
    private XYPlot mDayXYPlot = null;
    private SimpleXYSeries mAccelDataSeries = null;
    private Redrawer mRedrawer;

    double[] resultantBuffer = new double[1000]; //Approximately 1 minute of data. (16.666 * 60)
    int resultantBufferIndex = 0;
    int XValsIndex = 0;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_data_control);
        TabHost tabHost = (TabHost)findViewById(R.id.tabhost);
        tabHost.setup();
        TabHost.TabSpec daysTab = tabHost.newTabSpec("Day");
        daysTab.setContent(R.id.day);
        daysTab.setIndicator("Day");
        tabHost.addTab(daysTab);
        TabHost.TabSpec weeksTab = tabHost.newTabSpec("Week");
        weeksTab.setContent(R.id.week);
        weeksTab.setIndicator("Week");
        tabHost.addTab(weeksTab);
        TabHost.TabSpec monthTab = tabHost.newTabSpec("Month");
        monthTab.setContent(R.id.month);
        monthTab.setIndicator("Month");
        tabHost.addTab(monthTab);
        //Initialize Graph:
        mDayXYPlot = (XYPlot) findViewById(R.id.plotAccelDataDay);
        mDayXYPlot.setRangeStepMode(XYStepMode.INCREMENT_BY_VAL);
        mDayXYPlot.setDomainStepMode(XYStepMode.INCREMENT_BY_VAL);
        mDayXYPlot.setRangeBoundaries(0, 2, BoundaryMode.FIXED);
        mDayXYPlot.setRangeStepValue(2.00/5.0);
        mDayXYPlot.setDomainStepValue(DATA_HISTORY/5.0);
        mDayXYPlot.setDomainBoundaries(0,DATA_HISTORY,BoundaryMode.FIXED);
        mDayXYPlot.setDomainLabel("Time (mins)");
        mDayXYPlot.getDomainLabelWidget().pack();
        mDayXYPlot.setRangeLabel("Relative Activity Intensity");
        mDayXYPlot.getRangeLabelWidget().pack();
        mDayXYPlot.setRangeValueFormat(new DecimalFormat("#.###"));
        mDayXYPlot.setDomainValueFormat(new DecimalFormat("####.##"));
        mAccelDataSeries = new SimpleXYSeries("Activity Data");
        mAccelDataSeries.useImplicitXVals();
        LineAndPointFormatter accLine = new LineAndPointFormatter(Color.parseColor("#b30000"), null, null, null);
        mDayXYPlot.getDomainLabelWidget().getLabelPaint().setColor(Color.BLACK);
        mDayXYPlot.getDomainLabelWidget().getLabelPaint().setTextSize(32);
        mDayXYPlot.getRangeLabelWidget().getLabelPaint().setColor(Color.BLACK);
        mDayXYPlot.getRangeLabelWidget().getLabelPaint().setTextSize(32);
        mDayXYPlot.getGraphWidget().getDomainTickLabelPaint().setColor(Color.BLACK);
        mDayXYPlot.getGraphWidget().getRangeTickLabelPaint().setColor(Color.BLACK);
        mDayXYPlot.getGraphWidget().getDomainTickLabelPaint().setTextSize(36);
        mDayXYPlot.getGraphWidget().getRangeTickLabelPaint().setTextSize(36);
        mDayXYPlot.getGraphWidget().getDomainGridLinePaint().setColor(Color.WHITE);
        mDayXYPlot.getGraphWidget().getRangeGridLinePaint().setColor(Color.WHITE);
        mDayXYPlot.getLegendWidget().getTextPaint().setColor(Color.BLACK);
        mDayXYPlot.getLegendWidget().getTextPaint().setTextSize(32);
        mDayXYPlot.getTitleWidget().getLabelPaint().setTextSize(32);
        mDayXYPlot.getTitleWidget().getLabelPaint().setColor(Color.BLACK);
        mDayXYPlot.addSeries(mAccelDataSeries, accLine);
        Display display = ((WindowManager)getSystemService(Context.WINDOW_SERVICE)).getDefaultDisplay();
        final int refreshRate = (int)display.getRefreshRate();
        mRedrawer = new Redrawer(Arrays.asList(new Plot[]{mDayXYPlot}), refreshRate ,false);
    }

    @Override
    protected void onResume() {
        super.onResume();
        mRedrawer.start();
        try {
            readAllData();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            Log.e(TAG, "File Not Found!");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Override
    protected void onPause() {
        mRedrawer.pause();
        super.onPause();
    }

    @Override
    protected void onDestroy() {
        mRedrawer.finish();
        super.onDestroy();
    }

    private void readAllData() throws IOException {
            //Stop Writing Data in MainActivity.class:
        //1 search for CSVs that follow the format, and match the current day.
        File root = Environment.getExternalStorageDirectory();
        File directory = new File(root.getAbsolutePath() + "/MotionData");
        String fileStamp = "MotionSensorData_"+MainActivity.getDateStamp();
        File currentRawDataFile = new File(directory,fileStamp+".csv");
        //TODO: Search for file of classified data to analyze.
        File currentClassDataFile = new File(directory,fileStamp+"Classified.csv");
        if(currentRawDataFile.exists()) {
            //Fetch today's data
            csvReaderRawData = new CSVReader(new FileReader(currentRawDataFile));
            String[] nextLine;
            while ( (nextLine = csvReaderRawData.readNext()) !=null ) {
                if(nextLine.length>2) {
                    boolean[] isDouble = new boolean[3];
                    for (int i = 0; i < 3; i++) {
                        String tempString = nextLine[i];
                        isDouble[i] = (1==(tempString.length()-tempString.replace(".","").length())); //make sure there are no extra decimal places
                    }
                    if(isDouble[0] && isDouble[1] && isDouble[2]) {
                        double[] AccXYZ = {Double.parseDouble(nextLine[0]),
                                Double.parseDouble(nextLine[1]),
                                Double.parseDouble(nextLine[2])};
                        if(resultantBufferIndex<1000) {
                            resultantBuffer[resultantBufferIndex] = MainActivity.jniResultantAcc(AccXYZ);
                            resultantBufferIndex++;
                            adjustGraph();
                        } else {
                            //Get mean value and add to series.
                            double meanVal = MainActivity.jniMean(resultantBuffer, resultantBufferIndex);
                            plot(meanVal);
                            Log.d(TAG,"meanVal: "+String.valueOf(meanVal));
                            XValsIndex++;
                            resultantBufferIndex=0;
                        }
                    }
                }
            }
            //Calculate final minute of dataxz based on index.
            double meanVal = MainActivity.jniMean(resultantBuffer, resultantBufferIndex);
            //Plot data.
            plot(meanVal);
            adjustGraph();
        }
        //TODO: Read from data file with classified values.
        /*if(currentClassDataFile.exists()) {
            //Fetch class data.
        }*/
//        TODO: Retrieve & Queue Old Data for Viewing.
//        List<File> csvFiles = getCSVFiles(directory);
//        Log.e(TAG,"csvFilesSize: "+csvFiles.size());
//        Log.e(TAG,"csvFiles: "+csvFiles.toString());
    }

    /**
     * This 'plot' method adds a single data point to the series to be plotted in
     * the activity.
     * @param y Data value to be plotted. X vals are assigned implicitly.
     */
    private void plot(final double y) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                if(mAccelDataSeries.size()>DATA_HISTORY) {
                    mAccelDataSeries.removeFirst();
                }
                mAccelDataSeries.addLast(null,y);
            }
        });
    }

    private void adjustGraph() {
        double max = findGraphMax(mAccelDataSeries);
        double min = findGraphMin(mAccelDataSeries);
        if (max-min!=0) {
            mDayXYPlot.setRangeBoundaries(min, max, BoundaryMode.FIXED);
            mDayXYPlot.setRangeStepValue( (max-min) / 5.0 );
        } else {
            mDayXYPlot.setRangeBoundaries(min-0.004, max+0.004,BoundaryMode.AUTO);
            mDayXYPlot.setRangeStepValue((0.008+max-min)/5);
        }
        Number newMinX = Math.floor(0.00);
        Number newMaxX = Math.floor(XValsIndex);
        mDayXYPlot.setDomainBoundaries(newMinX, newMaxX, BoundaryMode.AUTO);
        mDayXYPlot.setDomainStepValue(XValsIndex/5.0);
    }

    private double findGraphMax(SimpleXYSeries s) {
        if (s.size() > 0 && s.getY(0)!=null) {
            double max = (double)s.getY(0);
            for (int i = 1; i < s.size(); i++) {
                double a = (double)s.getY(i);
                if(a>max) {
                    max = a;
                }
            }
            return max;
        } else
            return 0.0;
    }

    private double findGraphMin(SimpleXYSeries s) {
        if (s.size()>0 && s.getY(0)!=null) {
            double min = (double)s.getY(0);
            for (int i = 1; i < s.size(); i++) {
                double a = (double)s.getY(i);
                if(a<min) {
                    min = a;
                }
            }
            return min;
        } else {
            return 0.0;
        }
    }

    private List<File> getCSVFiles(File parentDirectory) {
        List<File> csvFiles = new ArrayList<>();
        Queue<File> files  = new LinkedList<>();
        files.addAll(Arrays.asList(parentDirectory.listFiles()));
        while (!files.isEmpty()) {
            File file = files.remove();
            if(file.isDirectory()) {
                files.addAll(Arrays.asList(file.listFiles()));
            } else if (file.getName().endsWith(".csv")) {
                csvFiles.add(file);
            }
        }
        return csvFiles;
    }

    private List<File> getMP3Files(File parentDirectory) {
        List<File> mp3Files = new ArrayList<>();
        Queue<File> files  = new LinkedList<>();
        files.addAll(Arrays.asList(parentDirectory.listFiles()));
        while (!files.isEmpty()) {
            File file = files.remove();
            if(file.isDirectory()) {
                files.addAll(Arrays.asList(file.listFiles()));
            } else if (file.getName().endsWith(".mp3")) {
                mp3Files.add(file);
            }
        }
        return mp3Files;
    }
}

package vcucmsc355.fitnessapp;

import java.util.ArrayList;
import android.content.Intent;
import android.app.IntentService;
import android.support.v4.content.LocalBroadcastManager;
import android.widget.ImageView;
import com.google.android.gms.location.ActivityRecognitionResult;
import com.google.android.gms.location.DetectedActivity;

/**
 * Created by Tony on 3/25/2017.
 */

public class ActivityRecognizedService extends IntentService {
    public static final String PACKAGE_NAME = "vcucmsc355.fitnessapp;";
    public static final String STRING_ACTION = PACKAGE_NAME + ".STRING_ACTION";
    public static final String STRING_EXTRA = PACKAGE_NAME + ".STRING_EXTRA";
    private ImageView image;
    public ActivityRecognizedService() {
        super("ActivityRecognize");
    }
    public ActivityRecognizedService(String name) {
        super(name);
    }
    @Override
    protected void onHandleIntent(Intent intent) {
        if(ActivityRecognitionResult.hasResult(intent)) {
            ActivityRecognitionResult result = ActivityRecognitionResult.extractResult(intent);
            Intent newIntent = new Intent(STRING_ACTION);
            ArrayList<DetectedActivity> detectedActivities = (ArrayList) result.getProbableActivities();
            newIntent.putExtra(STRING_EXTRA, detectedActivities);
            LocalBroadcastManager.getInstance(this).sendBroadcast(newIntent);
        }
    }
}

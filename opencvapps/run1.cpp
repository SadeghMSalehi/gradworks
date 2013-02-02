/**
 * opencv app
 *
 */

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    CvCapture *capture = cvCaptureFromAVI("/Users/joohwi/Dropbox/IntelInternship/myface.avi");
    if(!capture)
    {
        printf("!!! cvCaptureFromAVI failed (file not found?)\n");
        return -1;
    }

    int fps = (int) cvGetCaptureProperty(capture, CV_CAP_PROP_FPS);
    printf("* FPS: %d\n", fps);

    cvNamedWindow("display_video", CV_WINDOW_AUTOSIZE);

    IplImage* frame = NULL;
    char key = 0;

    while (key != 'q')
    {
        frame = cvQueryFrame(capture);
        if (!frame)
        {
            printf("!!! cvQueryFrame failed: no frame\n");
            break;
        }

        cvShowImage("display_video", frame);
        
        key = cvWaitKey(1000 / fps);
    }
    
    cvReleaseCapture(&capture);
    cvDestroyWindow("display_video");
}
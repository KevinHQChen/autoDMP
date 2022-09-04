#include "improc/improc.hpp"

#if 0
void imgAnnotation(QueueFPS<cv::Mat>& imgQueue, cv::Mat& img) {
            // adding basic frame annotation
            std::ostringstream label;
            label << std::fixed << std::setprecision(2) << imgQueue.getFPS();

            std::ostringstream label2;
            label2 << imgQueue.counter_();

            cv::putText(img, label.str(), cv::Point(0, 60), cv::FONT_HERSHEY_SIMPLEX, 0.25, cv::Scalar::all(0));
            cv::putText(img, label2.str(), cv::Point(0, 75), cv::FONT_HERSHEY_SIMPLEX, 0.25, cv::Scalar::all(0));
}
#endif

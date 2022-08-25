#pragma once

#include "util/util.hpp"
#include "cam/cam.hpp"


#define NUM_TEMPLATES 4

struct pose {
    int rotation;
    cv::Point matchLoc;
};

struct ChannelPose {
    std::vector<int> rotAngle;
    std::vector<cv::Rect> chanBBox;
    std::vector<cv::Rect> rotChanBBox;
};

void matInfo(cv::Mat mat);

void imgAnnotation(QueueFPS<cv::Mat>& imgQueue, cv::Mat& img);

void rotateMat(cv::Mat& src, cv::Mat& dst, double angle);

void rotateMatCropped(cv::Mat& src, cv::Mat& dst, double angle);

int tmSetup(config conf,
            cam* onlineCam,
            cv::VideoCapture* offlineCam,
            std::array<cv::Mat, NUM_TEMPLATES>& templateImg,
            std::array<cv::cuda::GpuMat, NUM_TEMPLATES>& templateImgGPU,
            ChannelPose& chanPose);

void imCap(config conf,
           QueueFPS<cv::Mat>& rawFramesQueue,
           QueueFPS<cv::Mat>& preFramesQueue,
           cam* onlineCam,
           cv::VideoCapture* offlineCam,
           bool& run);

void imProc(config conf,
            std::array<cv::Mat, NUM_TEMPLATES>& templateImg,
            std::array<cv::cuda::GpuMat, NUM_TEMPLATES>& templateImgGPU,
            ChannelPose& chanPose,
            QueueFPS<cv::Mat>& preFramesQueue,
            std::vector<QueueFPS<cv::Mat>*>& tempResultsQueues,
            std::vector<QueueFPS<cv::Mat>*>& procFramesQueues,
            std::vector<QueueFPS<cv::Point>*>& procDataQueues,
            bool& run);

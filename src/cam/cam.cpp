#include "cam/cam.hpp"

using namespace std;
using namespace spdlog;


cam::cam(int cameraIdx) : cameraIndex(cameraIdx) {
    // init libs
    returnCode = AT_InitialiseLibrary();
    if (returnCode != AT_SUCCESS) {
        error("FAIL: zyla can not initialise library, return {}\n", returnCode);
        return;
    } else
        info("zyla initialized library\n");
    returnCode = AT_InitialiseUtilityLibrary();
    if (returnCode != AT_SUCCESS) {
        error("FAIL: zyla can not initialise utility library, return {}\n", returnCode);
        return;
    } else
        info("zyla initialized utility library\n");

    // open handle
    returnCode = AT_Open(cameraIndex, &handle);
    if (returnCode != AT_SUCCESS) {
        error("FAIL: zyla can not open, return {}\n", returnCode);
        return;
    } else
        info("zyla opened with handle {}\n", handle);
}

cam::~cam() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    // close handle
    returnCode = AT_Close(handle);
    if (returnCode != AT_SUCCESS) {
        mOut << "FAIL: zyla can not close, return " <<
            returnCode << "\n";
    } else {
        mOut << "zyla closed" << "\n";
    }

    // clean up libs
    returnCode = AT_FinaliseLibrary();
    if (returnCode != AT_SUCCESS) {
        mOut << "FAIL: zyla can not finalise library, return " <<
            returnCode << "\n";
    } else {
        mOut << "zyla finalized library" << "\n";
    }
    returnCode = AT_FinaliseUtilityLibrary();
    if (returnCode != AT_SUCCESS)
    {
        mOut << "FAIL: zyla can not finalise utility library, return " <<
            returnCode << "\n";
    }
    else
    {
        mOut << "zyla finalized utility library" << "\n";
    }
}

void cam::get(camSettings &s) {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    //// bool
    for (s.boolMapIterator = s.boolMap.begin();
        s.boolMapIterator != s.boolMap.end();
        s.boolMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.boolMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.boolMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.boolMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.boolMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not get " << s.boolMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.boolMapIterator->first;
            if (readable)
            {
                AT_BOOL value;
                returnCode = AT_GetBool(handle, s.boolMapIterator->first,
                    &value);
                s.boolMapIterator->second = value;
                mOut << " = " << value;
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
    //// int
    for (s.intMapIterator = s.intMap.begin();
        s.intMapIterator != s.intMap.end();
        s.intMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.intMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.intMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.intMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.intMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not get " << s.intMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.intMapIterator->first;
            if (readable)
            {
                AT_64 value;
                returnCode = AT_GetInt(handle, s.intMapIterator->first,
                    &value);
                s.intMapIterator->second = value;
                mOut << " = " << value;
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
    //// float
    for (s.floatMapIterator = s.floatMap.begin();
        s.floatMapIterator != s.floatMap.end();
        s.floatMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.floatMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.floatMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.floatMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.floatMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not get " << s.floatMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.floatMapIterator->first;
            if (readable)
            {
                double value;
                returnCode = AT_GetFloat(handle, s.floatMapIterator->first,
                    &value);
                s.floatMapIterator->second = value;
                mOut << " = " << value;
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
    //// enum
    for (s.enumMapIterator = s.enumMap.begin();
        s.enumMapIterator != s.enumMap.end();
        s.enumMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.enumMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.enumMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.enumMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.enumMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not get " << s.enumMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.enumMapIterator->first;
            if (readable)
            {
                int index = 0;
                AT_WC *str = new AT_WC[256];
                returnCode = AT_GetEnumIndex(handle,
                    s.enumMapIterator->first,
                    &index);
                returnCode = AT_GetEnumStringByIndex(handle,
                    s.enumMapIterator->first,
                    index, str, 256);
                //delete[] s.enumMapIterator->second; // small leak
                s.enumMapIterator->second = str;
                mOut << " = " << str;
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
}

void cam::set(camSettings &s) {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    //// bool
    for (s.boolMapIterator = s.boolMap.begin();
        s.boolMapIterator != s.boolMap.end();
        s.boolMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.boolMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.boolMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.boolMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.boolMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not set " << s.boolMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.boolMapIterator->first;
            if (writable)
            {
                returnCode = AT_SetBool(handle,
                    s.boolMapIterator->first,
                    s.boolMapIterator->second);
                if (returnCode != AT_SUCCESS)
                    mOut << " can NOT set";
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
    //// int
    for (s.intMapIterator = s.intMap.begin();
        s.intMapIterator != s.intMap.end();
        s.intMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.intMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.intMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.intMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.intMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not set " << s.intMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.intMapIterator->first;
            if (writable)
            {
                returnCode = AT_SetInt(handle,
                    s.intMapIterator->first,
                    s.intMapIterator->second);
                if (returnCode != AT_SUCCESS)
                    mOut << " can NOT set";
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
    //// float
    for (s.floatMapIterator = s.floatMap.begin();
        s.floatMapIterator != s.floatMap.end();
        s.floatMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.floatMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.floatMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.floatMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.floatMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not set " << s.floatMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.floatMapIterator->first;
            if (writable)
            {
                returnCode = AT_SetFloat(handle,
                    s.floatMapIterator->first,
                    s.floatMapIterator->second);
                if (returnCode != AT_SUCCESS)
                    mOut << " can NOT set";
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
    //// enum
    for (s.enumMapIterator = s.enumMap.begin();
        s.enumMapIterator != s.enumMap.end();
        s.enumMapIterator++)
    {
        AT_BOOL implemented;
        returnCode = AT_IsImplemented(handle, s.enumMapIterator->first,
            &implemented);
        AT_BOOL readOnly;
        returnCode = AT_IsReadOnly(handle, s.enumMapIterator->first,
            &readOnly);
        AT_BOOL readable;
        returnCode = AT_IsReadable(handle, s.enumMapIterator->first,
            &readable);
        AT_BOOL writable;
        returnCode = AT_IsWritable(handle, s.enumMapIterator->first,
            &writable);
        if (returnCode != AT_SUCCESS)
        {
            mOut << "FAIL: can not set " << s.enumMapIterator->first << "\n";
        }
        else
        {
            mOut << "i" << implemented << "o" << readOnly <<
                "r" << readable << "w" << writable << "\t" <<
                s.enumMapIterator->first;
            if (writable)
            {
                returnCode = AT_SetEnumString(handle,
                    s.enumMapIterator->first,
                    s.enumMapIterator->second);
                if (returnCode != AT_SUCCESS)
                    mOut << " can NOT set";
            }
            mOut << "\n";
        }
    }
    mOut << "\n";
}

void cam::start(const int &Ts) {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    //// flush queue and wait buffers;
    returnCode = AT_Flush(handle);

    //// get info
    samplePeriod = Ts;
    // get size of memory to be allocated for storing each acquired image
    returnCode = AT_GetInt(handle, L"ImageSizeBytes", &imageSizeBytes);
    // cast from AT_64 to int so this value can be used in AT_QueueBuffer
    bufferSize = (int)(imageSizeBytes);

    returnCode = AT_GetFloat(handle, L"FrameRate", &frameRate);
    // set the circular buffer size (number of buffers)
    queueLength = (int)(frameRate * samplePeriod + 1); // min length = 1;
    returnCode = AT_GetInt(handle, L"AOIStride", &imageStride);
    returnCode = AT_GetInt(handle, L"AOIWidth", &imageWidth);
    returnCode = AT_GetInt(handle, L"AOILeft", &imageLeft);
    returnCode = AT_GetInt(handle, L"AOIHeight", &imageHeight);
    returnCode = AT_GetInt(handle, L"AOITop", &imageTop);
    mOut << "Frame Size (bytes): " << imageSizeBytes << "\n"
            << "Framerate: " << frameRate << "\n"
            << "ROI width (pixels): " << imageWidth << "\n"
            << "ROI left (pixels): " << imageLeft << "\n"
            << "ROI height (pixels): " << imageHeight << "\n"
            << "ROI top (pixels): " << imageTop << "\n"
            << "ROI Stride (bytes): " << imageStride << "\n"
            << "Queue Length: " << queueLength << "\n";

    int index = 0;    // [0,1,2,3]:[Mono16,Mono12,Mono12Packed,Mono32]
    // AT_WC = andor-specific 16-bit wide character type
    // (used to represent all feature names, enumerated options, and string feature values)
    AT_WC *str = new AT_WC[256];
    // get pixel encoding value index [0,1,2,3]:[Mono16,Mono12,Mono12Packed,Mono32]
    returnCode = AT_GetEnumIndex(handle, L"Pixel Encoding", &index);
    // convert index to string representation
    returnCode = AT_GetEnumStringByIndex(handle, L"Pixel Encoding", index, str, 256);
    imageEncode = str;

    //// allocate buffer
    buffers = new unsigned char*[queueLength];
    alignedBuffers = new unsigned char*[queueLength];
    // memory allocated should be aligned on an 8-byte (64-bit) boundary
    // this helps system performance and prevents alignment faults
    for (int i = 0; i < queueLength; i++) {
        // each index of buffers*[] is a uc pointer
        // so for each index, we'll point it to a newly allocated space
        // that holds a single image divided into uc-sized chunks
        buffers[i] = new unsigned char[bufferSize + 8]; // add 7 to allow data alignment
        // this magic somehow aligns it to 64-bit
        alignedBuffers[i] = reinterpret_cast<unsigned char*>
            ((reinterpret_cast<unsigned    long long>(buffers[i% queueLength]) + 7) & ~7);
    }
    //// pass buffers to fifo queue
    // AT_QueueBuffer will configure the corresponding area of memory to store the acquired images
    for (int i = 0; i < queueLength; i++) {
        returnCode = AT_QueueBuffer(handle, alignedBuffers[i], bufferSize);
    }
    //// start acquisition
    returnCode = AT_Command(handle, L"AcquisitionStart");
    std::cerr << "acquisition start returns " << returnCode << "\n";
    accumNumFrames = 0;
}

void cam::stop() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    //// stop acquisition
    returnCode = AT_Command(handle, L"AcquisitionStop");
    std::cerr << "acquisition stop returns " << returnCode << "\n";
    returnCode = AT_Flush(handle);
    std::cerr << "queue/wait buffer flush returns " << returnCode << "\n";
    //// free buffer
    for (int i = 0; i < queueLength; i++) {
        delete[] buffers[i];
    }
    delete[] buffers;
    delete[] alignedBuffers;
}

int cam::process(cv::Mat &image) {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    // grab buffer
    unsigned char* pointer;
    int size;
    // if a frame is available, store its address in pointer, size in size, and timeout in ms
    // AT_WaitBuffer will put the calling thread to sleep until timeout elapses or an image becomes available
    // (if timeout equals 0 it will only get existing frames without waiting)
    // (0 timeout causes AT_WaitBuffer to hang quite often, setting to 15ms helps a lot)
    // std::cerr << "zyla process wait buffer returns ";
    returnCode = AT_WaitBuffer(handle, &pointer, &size, 15);
    // std::cerr << returnCode << "\n";

    if (returnCode == 0) {
        // re-queue circular buffer
        returnCode = AT_QueueBuffer(handle, alignedBuffers[accumNumFrames % queueLength], bufferSize);
        // std::cerr << "re-queue returns " << returnCode << "\n";
        accumNumFrames++;
        std::cerr << "accumNumFrames: " << accumNumFrames << "\n";
        // clean up buffer
        image = cv::Mat(imageHeight, imageWidth, CV_16UC1);
        returnCode = AT_ConvertBuffer(pointer, reinterpret_cast<unsigned char*>(image.data),
            imageWidth, imageHeight, imageStride, imageEncode, L"Mono16");
        // std::cerr << "convert returns " << returnCode << "\n";
    }
    return returnCode;
}

#pragma once

#include <iostream>
#include <opencv2\opencv.hpp>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <windows.h>
#include <sstream>
#include <string>

using namespace std;
using namespace cv;

#define N 200
#define InvDim 64
#define GHMIDir "F:\\Gaussian-Hermite3D\\FaceGHMIs\\FaceGHMIs\\"
#define IMGDir "F:\\DatabaseFace\\CASIA-3D-FaceV1\\CASIA-3D-FaceV1\\3D-Face-BMP\\"

typedef double Invtype;

class Sort
{
public:
	int m_imageCount;
public:
	Sort(String fname, double ChiDist[N][N], vector<string> &names);
	~Sort();

	/**********Sort image by their Dist value, and display them.*********/
	void Sort::SortAndDisplay(double Dist[N][N], int imagenumber, vector<string> &names);
private:
	void MultiImage_OneWin(const std::string& MultiShow_WinName, const vector<Mat>& SrcImg_V, Size SubPlot, Size ImgMax_Size);
	string Trim(string& str);
	int Find(double targetvalue, double OriginArray[]);
	void ReadDotCsv(String fname, double ChiDist[N][N], vector<string> &names);
	void ReadOrderList(String filename, vector<string> &names, int nums[N]);
};
#include "sort.h"

Sort::Sort(String fname, double Dist[N][N], vector<string> &names)
{
	this->m_imageCount = 0;
	ReadDotCsv(fname, Dist, names);
}

Sort::~Sort()
{

}

void Sort::MultiImage_OneWin(const string& MultiShow_WinName, const vector<Mat>& SrcImg_V, Size SubPlot, Size ImgMax_Size)
{
	//Reference : http://blog.csdn.net/yangyangyang20092010/article/details/21740373  

	//************* Usage *************//  
	//vector<Mat> imgs(4);  
	//imgs[0] = imread("F:\\SA2014.jpg");  
	//imgs[1] = imread("F:\\SA2014.jpg");  
	//imgs[2] = imread("F:\\SA2014.jpg");  
	//imgs[3] = imread("F:\\SA2014.jpg");  
	//MultiImage_OneWin("T", imgs, cvSize(2, 2), cvSize(400, 280));  

	//Window's image  
	Mat Disp_Img;
	//Width of source image  
	Size Img_OrigSize(SrcImg_V[0].cols, SrcImg_V[0].rows);
	//******************** Set the width for displayed image ********************//  
	//Width vs height ratio of source image  
	float WH_Ratio_Orig = Img_OrigSize.width / (float)Img_OrigSize.height;
	Size ImgDisp_Size(100, 100);
	if (Img_OrigSize.width > ImgMax_Size.width)
	{
		ImgDisp_Size.width = ImgMax_Size.width;
		ImgDisp_Size.height = (int)ImgMax_Size.width / WH_Ratio_Orig;
	}
	else if (Img_OrigSize.height > ImgMax_Size.height)
	{
		ImgDisp_Size.width = (int)ImgMax_Size.height*WH_Ratio_Orig;
		ImgDisp_Size.height = ImgMax_Size.height;
	}
	else
	{
		ImgDisp_Size.width = Img_OrigSize.width;
		ImgDisp_Size.height = Img_OrigSize.height;
	}
	//******************** Check Image numbers with Subplot layout ********************//  
	int Img_Num = (int)SrcImg_V.size();
	if (Img_Num > SubPlot.width * SubPlot.height)
	{
		cout << "Your SubPlot Setting is too small !" << endl;
		exit(0);
	}
	//******************** Blank setting ********************//  
	Size DispBlank_Edge(80, 60);
	Size DispBlank_Gap(15, 15);
	//******************** Size for Window ********************//  
	Disp_Img.create(Size(ImgDisp_Size.width*SubPlot.width + DispBlank_Edge.width + (SubPlot.width - 1)*DispBlank_Gap.width,
		ImgDisp_Size.height*SubPlot.height + DispBlank_Edge.height + (SubPlot.height - 1)*DispBlank_Gap.height), CV_8UC3);
	Disp_Img.setTo(0);//Background  
	//Left top position for each image  
	int EdgeBlank_X = (Disp_Img.cols - (ImgDisp_Size.width*SubPlot.width + (SubPlot.width - 1)*DispBlank_Gap.width)) / 2;
	int EdgeBlank_Y = (Disp_Img.rows - (ImgDisp_Size.height*SubPlot.height + (SubPlot.height - 1)*DispBlank_Gap.height)) / 2;
	Point LT_BasePos(EdgeBlank_X, EdgeBlank_Y);
	Point LT_Pos = LT_BasePos;

	//Display all images  
	for (int i = 0; i < Img_Num; i++)
	{
		//Obtain the left top position  
		if ((i%SubPlot.width == 0) && (LT_Pos.x != LT_BasePos.x))
		{
			LT_Pos.x = LT_BasePos.x;
			LT_Pos.y += (DispBlank_Gap.height + ImgDisp_Size.height);
		}
		//Writting each to Window's Image  
		Mat imgROI = Disp_Img(Rect(LT_Pos.x, LT_Pos.y, ImgDisp_Size.width, ImgDisp_Size.height));
		resize(SrcImg_V[i], imgROI, Size(ImgDisp_Size.width, ImgDisp_Size.height));

		LT_Pos.x += (DispBlank_Gap.width + ImgDisp_Size.width);
	}

	//Get the screen size of computer  
	int Scree_W = GetSystemMetrics(SM_CXSCREEN);
	int Scree_H = GetSystemMetrics(SM_CYSCREEN);

	namedWindow(MultiShow_WinName.c_str(), WINDOW_AUTOSIZE);
	moveWindow(MultiShow_WinName.c_str(), (Scree_W - Disp_Img.cols) / 2, (Scree_H - Disp_Img.rows) / 2);//Centralize the window
	imshow(MultiShow_WinName.c_str(), Disp_Img);
	//waitKey(0);
}

string Sort::Trim(string& str)
{
	str.erase(0, str.find_first_not_of(" \t\r\n"));
	str.erase(str.find_last_not_of(" \t\r\n") + 1);
	return str;
}

//此处可根据算法改进
int Sort::Find(double targetvalue, double OriginArray[])
{
	for (int k = 0; k < this->m_imageCount; k++)
	if (targetvalue == OriginArray[k]) return k + 1;
}

//此处.csv是一个下三角矩阵，输出Dist是一个imageCount x imageCount的对称矩阵。
void Sort::ReadDotCsv(String fname, double Dist[N][N], vector<string> &names)
{
	const char* filename = fname.c_str();
	ifstream fin(filename);
	string s;

	string line;
	//getline(fin, line);

	//int index;

	//index = line.find_last_of(',');
	//// firstname = "Mark" lastname = "Tompkin"  
	//this->m_imageCount = atof(line.substr(index + 1).c_str());

	int i = 1;
	vector< vector<Invtype> > Invs;
	vector<Invtype> Inv;

	while (getline(fin, line)) {
		//cout << line << endl;

		istringstream sin(line);
		vector<string> fields;
		string field;
		/*getline(sin, field, ',');
		int i = atof(field.c_str()) - 1;*/

		while (getline(sin, field, ','))
		{
			fields.push_back(field);
		}

		names.push_back((string)Trim(fields[0]).c_str());
		for (int j = 0; j < InvDim; j++)
		{
			Inv.push_back(atof(Trim(fields[j+1]).c_str()));
		}
		Invs.push_back(Inv);
		Inv.clear();
		i++;
		if (i == 201) 
			break;
	}

	this->m_imageCount = i-1;

	for (size_t filenum1 = 0; filenum1 != Invs.size(); filenum1++)
	{
		for (size_t filenum0 = 0; filenum0 < filenum1; filenum0++)
		{
			vector<Invtype> Inv0, Inv1;
			Inv0 = Invs[filenum0];
			Inv1 = Invs[filenum1];

			for (int i = 0; i < InvDim; i++)
			{
				if (0 == abs(Inv0[i] + Inv1[i]))
					Dist[filenum0][filenum1] = (Inv0[i] - Inv1[i])*(Inv0[i] - Inv1[i]);
				else
					Dist[filenum0][filenum1] += (Inv0[i] - Inv1[i])*(Inv0[i] - Inv1[i]) / abs(Inv0[i] + Inv1[i]);

				Dist[filenum1][filenum0] = Dist[filenum0][filenum1];
			}
		}
	}
}

void Sort::ReadOrderList(String filename, vector<string> &names, int nums[N])
{
	string newline;
	//int currentfile;

	ifstream pathOrder(filename.c_str());

	string name;
	string num;

	int idx = 0;
	while (getline(pathOrder, newline))
	{
		int index;

		index = newline.find_last_of('\t');
		// firstname = "Mark" lastname = "Tompkin"  
		name = newline.substr(0, index);
		num = newline.substr(index + 1);

		nums[idx] = atof(num.c_str());

		names.push_back(name);
		idx++;
	}
}

void Sort::SortAndDisplay(double Dist[N][N], int imagenumber, vector<string> &names)
{
	double a[N], temp[N];
	for (int n = 0; n < this->m_imageCount; n++)
	{
		a[n] = Dist[imagenumber - 1][n];
		temp[n] = a[n];
	}

	sort(a, a + this->m_imageCount);

	int order[N];
	for (int n = 0; n < this->m_imageCount; n++)
	{
		order[n] = Find(a[n], temp);
		cout << order[n] << "\t" << endl;
	}

	//String orderstr = GHMIDir + (String)"order.txt";

	int nums[N];

	/*ReadOrderList(orderstr, names, nums);*/

	vector<Mat> imgs(30);

	
	for (int i = 0; i < 30; i++)
	{
		int index = 0;//这里要对index赋0不然出错。
		String dirwrl = names[order[i] - 1];
		
		index = dirwrl.find_first_of('L',index);

		String folder = dirwrl.substr(index+2, 3);

		index = dirwrl.find_first_of('.', index);
		String file = dirwrl.substr(index - 7, 7);
		// firstname = "Mark" lastname = "Tompkin"  
		

		String dirbmp;
		dirbmp = IMGDir + folder + "\\" + file + ".bmp";
		imgs[i] = imread(dirbmp);
	}

	MultiImage_OneWin("Sort Results", imgs, Size(6, 5), Size(100, 100));
}
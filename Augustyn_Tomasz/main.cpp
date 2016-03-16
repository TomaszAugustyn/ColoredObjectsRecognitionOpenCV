#include "opencv2/opencv.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>


using namespace cv;
using namespace std;

#define M_PI           3.14159265358979323846  /* pi */

struct AngleAndElemNumber
{
	double dPolarAngle;
	int iElementNumber;
};

Point2f CalculateCenterPointOfPoints( std::vector<Point2f> mc2 )
{
	size_t iNumberOfPoints = mc2.size();
	float dCenterX = 0, dCenterY = 0;
	Point2f oCenterPoint;

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		dCenterX = dCenterX + mc2[i].x;
		dCenterY = dCenterY + mc2[i].y;
	}
	dCenterX = dCenterX / iNumberOfPoints;
	dCenterY = dCenterY / iNumberOfPoints;
	oCenterPoint.x = dCenterX;
	oCenterPoint.y = dCenterY;

	return oCenterPoint;

}

std::vector<Point2f> RelateVectorOfPointsToTheCenterPoint( std::vector<Point2f> mc2, Point2f oCenterPoint )
{
	std::vector<Point2f> mc2_related;
	size_t iNumberOfPoints = mc2.size();
	Point2f temp_point;

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		temp_point.x = mc2[i].x - oCenterPoint.x;
		temp_point.y = mc2[i].y - oCenterPoint.y;
		mc2_related.push_back( temp_point );
	}
	return mc2_related;

}

std::vector<AngleAndElemNumber> ConvertToPolarCoordinates(std::vector<Point2f> mc2_related)
{
	int iNumberOfPoints = mc2_related.size();
	std::vector<AngleAndElemNumber> VectorOfAngleAndElem;
	AngleAndElemNumber TempStruct;
	float tempX = 0.0, tempY = 0.0;

	for ( int i = 0; i < iNumberOfPoints; i++ )
	{
		TempStruct.iElementNumber = i;
		tempX = mc2_related[i].x;
		tempY = mc2_related[i].y;

		if ( ( tempX > 0 ) && ( tempY >= 0 ) )
		{
			TempStruct.dPolarAngle = atan(tempY / tempX);
		}
		else if ( ( tempX > 0) && (tempY < 0 ) )
		{
			TempStruct.dPolarAngle = atan(tempY / tempX) + 2 * M_PI;
		}
		else if ( tempX < 0 )
		{
			TempStruct.dPolarAngle = atan(tempY / tempX) + M_PI;
		}
		else if ( (tempX == 0) && (tempY > 0) )
		{
			TempStruct.dPolarAngle = M_PI / 2;
		}
		else if ( (tempX == 0) && (tempY < 0) )
		{
			TempStruct.dPolarAngle = 3 * M_PI / 2;
		}
		else
		{
			TempStruct.dPolarAngle = 0;
		}

		VectorOfAngleAndElem.push_back( TempStruct );
	}

	return VectorOfAngleAndElem;
}

int Partition( std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r ) // dzielimy tablice na dwie czesci, w pierwszej wszystkie liczby sa mniejsze badz rowne x, w drugiej wieksze lub rowne od x
{
	double x = VectorOfAngleAndElem[p].dPolarAngle; // obieramy x
	int i = p, j = r, index; // i, j - indeksy w tabeli
	double w;
	while (true) // petla nieskonczona - wychodzimy z niej tylko przez return j
	{
		while (VectorOfAngleAndElem[j].dPolarAngle > x) // dopoki elementy sa wieksze od x
			j--;
		while (VectorOfAngleAndElem[i].dPolarAngle < x) // dopoki elementy sa mniejsze od x
			i++;
		if (i < j) // zamieniamy miejscami gdy i < j
		{
			w = VectorOfAngleAndElem[i].dPolarAngle;
			index = VectorOfAngleAndElem[i].iElementNumber;

			VectorOfAngleAndElem[i].dPolarAngle = VectorOfAngleAndElem[j].dPolarAngle;
			VectorOfAngleAndElem[i].iElementNumber = VectorOfAngleAndElem[j].iElementNumber;

			VectorOfAngleAndElem[j].dPolarAngle = w;
			VectorOfAngleAndElem[j].iElementNumber = index;
			i++;
			j--;
		}
		else // gdy i >= j zwracamy j jako punkt podzialu tablicy
			return j;
	}
}

void Quicksort( std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r ) // sortowanie szybkie
{
	int q;
	if (p < r)
	{
		q = Partition(VectorOfAngleAndElem, p, r); // dzielimy tablice na dwie czesci; q oznacza punkt podzialu
		Quicksort(VectorOfAngleAndElem, p, q); // wywolujemy rekurencyjnie quicksort dla pierwszej czesci tablicy
		Quicksort(VectorOfAngleAndElem, q + 1, r); // wywolujemy rekurencyjnie quicksort dla drugiej czesci tablicy
	}
}

std::vector<Point2f> EnumerateVerticiesClockwise( std::vector<AngleAndElemNumber> VectorOfAngleAndElem, std::vector<Point2f> mc2 )
{
	size_t iNumberOfPoints = VectorOfAngleAndElem.size();
	std::vector<Point2f> mc2_clockwise;
	Point2f TempPoint;
	int ElementNr;

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		ElementNr = VectorOfAngleAndElem[i].iElementNumber;
		TempPoint.x = mc2[ElementNr].x;
		TempPoint.y = mc2[ElementNr].y;
		mc2_clockwise.push_back( TempPoint );
	}

	return mc2_clockwise;
}

int main(int, char)
{
	//odczyt pliku tekstowego
	vector<string> nazwy_obrazkow;
	string pojedyncza_linia;
	fstream plik;
	plik.open("nazwy_zdjec/nazwy_zdjec.txt", ios::in);
	if ( plik.good() == false )
	{
		cout << "Podany plik nie istnieje lub podales bledna sciezke!" << endl;
		exit(0);
	}
	while ( getline( plik, pojedyncza_linia ) )
	{
		nazwy_obrazkow.push_back( pojedyncza_linia );
	}
	plik.close();


	//ladowanie obrazkow do wektora Mat i skalowanie ich
	vector <Mat> frames;
	Mat temp_frame;
	for ( size_t i = 0; i < nazwy_obrazkow.size(); i++ )
	{
		temp_frame = imread("zdjecia/" + nazwy_obrazkow[i] , CV_LOAD_IMAGE_GRAYSCALE);
		Mat resized;

		if (temp_frame.rows <= temp_frame.cols)
		{
			resized.create(1170, 2080, temp_frame.type());
			
		}
		else
		{
			Mat resized(2080, 1170, temp_frame.type());
		}
		
		cv::resize(temp_frame, resized, resized.size());
		frames.push_back(resized);
	}



	Mat frame, threshold_img, element, dilation_dst, dest_img, transform_matrix, erode_dst;
	int dilation_size = 3;
	vector< vector <Point> > contours;
	namedWindow("okno", 1);
	namedWindow("okno2", 1);
	frame = imread("zelki.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	threshold(frames[0], threshold_img, 30, 255, THRESH_BINARY);

	element = getStructuringElement(MORPH_ELLIPSE,
		Size(2 * dilation_size + 1, 2 * dilation_size + 1),
		Point(dilation_size, dilation_size));

	dilate(threshold_img, dilation_dst, element, cv::Point(-1,-1), 3);
	erode(dilation_dst, erode_dst, element);

	findContours(erode_dst, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);

	Scalar color(255, 255, 255);
	for (int i = 1; i < contours.size(); i++)
	{
		//drawContours(dilation_dst, contours, i, color, CV_FILLED);
		drawContours(erode_dst, contours, i, color, 2, 8);
	}

	/// Get the moments
	vector<Moments> mu(contours.size());
	for (int i = 0; i < contours.size(); i++)
	{
		mu[i] = moments(contours[i], false);
	}

	///  Get the mass centers:
	vector<Point2f> mc(contours.size());
	vector<Point2f> mc2;
	vector<Point2f> mc2_related;
	Point2f oCenterPoint;
	vector<AngleAndElemNumber> VectorOfAngleAndElem;
	vector<Point2f> mc2_clockwise;


	for (int i = 0; i < contours.size(); i++)
	{
		mc[i] = Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
	}

	for (int i = 1; i < mc.size(); i++)
	{
		mc2.push_back(mc[i]);
	}

	oCenterPoint = CalculateCenterPointOfPoints( mc2 );
	mc2_related = RelateVectorOfPointsToTheCenterPoint( mc2, oCenterPoint );
	VectorOfAngleAndElem = ConvertToPolarCoordinates( mc2_related );
	Quicksort( VectorOfAngleAndElem, 0, (int)VectorOfAngleAndElem.size() - 1 );
	mc2_clockwise = EnumerateVerticiesClockwise( VectorOfAngleAndElem, mc2 );

	vector<Point2f> punktyDocelowe;
	Point2f p1, p2, p3, p4;

	p4.x = 800;
	p4.y = 0;

	p3.x = 0;
	p3.y = 0;

	p2.x = 0;
	p2.y = 500;

	p1.x = 800;
	p1.y = 500;



	punktyDocelowe.push_back(p1);
	punktyDocelowe.push_back(p2);
	punktyDocelowe.push_back(p3);
	punktyDocelowe.push_back(p4);

	transform_matrix = getPerspectiveTransform(mc2_clockwise, punktyDocelowe);
	warpPerspective(frames[0], dest_img, transform_matrix, Size(800, 500));
	imshow("okno2", erode_dst);
	imshow("okno", dest_img);

	waitKey();
	return 0;
}
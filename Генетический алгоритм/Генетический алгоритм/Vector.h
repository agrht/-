#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <ctime>
const int MAX_VECTOR_SIZE = 100000000;

using namespace std;
template <class ValType>
class TVector
{
	/*struct 
	{
		Type x, y;
	};*/
protected:
	ValType *pVector;
	int Size;       // ������ �������
	int StartIndex; // ������ ������� �������� �������
public:
	TVector(int a = 10, int b = 0);
	TVector(const TVector &a);                // ����������� �����������
	~TVector();
	int GetSize() { return Size; } // ������ �������
	int GetStartIndex() { return StartIndex; } // ������ ������� ��������
	ValType& operator[](int pos);             // ������
	bool operator==(const TVector &a) const;  // ���������
	bool operator!=(const TVector &a) const;  // ���������
	TVector& operator=(const TVector &a);     // ������������
	// ��������� ��������
	TVector  operator+(const ValType &val);   // ��������� ������
	TVector  operator-(const ValType &val);   // ������� ������
	TVector  operator*(const ValType &val);   // �������� �� ������
	TVector  operator*(const int &val);   // �������� �� ������
	// ��������� ��������
	TVector  operator+(const TVector &v);     // ��������
	TVector  operator-(const TVector &v);     // ���������
	ValType  operator*(const TVector &v);     // ��������� ������������
	//Type operator %() const { return pow(x*x + y * y, 0.5); };
	// ����-�����
	friend istream& operator>>(istream &is, TVector &a)
	{
		for (int i = 0; i < a.Size; i++)
			is >> a.pVector[i];
		return is;
	}
	friend ostream& operator<<(ostream &os, const TVector &a)
	{
		for (int i = 0; i < a.Size; i++)
			os << a.pVector[i] << ' ';
		return os;
	}
};

template <class ValType> //�����������
TVector<ValType>::TVector(int a, int b)
{
	if (a < 0 || b < 0 || a > MAX_VECTOR_SIZE)
	{
		throw - 1;
	}
	Size = a;
	StartIndex = b;
	pVector = new ValType[Size];
	for (int i = 0; i < Size; i++)
	{
		pVector[i] = 0;
	}
}

template <class ValType> //����������� �����������
TVector<ValType>::TVector(const TVector<ValType> &a)
{
	Size = a.Size;
	StartIndex = a.StartIndex;
	pVector = new ValType[Size];
	for (int i = 0; i < Size; i++)
	{
		pVector[i] = a.pVector[i];
	}
}

template <class ValType>
TVector<ValType>::~TVector()
{
	delete[] pVector;
}

template <class ValType> // ������
ValType& TVector<ValType>::operator[](int pos)
{
	if (pos - StartIndex < 0 || pos - StartIndex >= Size)
		throw - 1;
	return pVector[pos - StartIndex];
}

template <class ValType> // ���������
bool TVector<ValType>::operator==(const TVector &a) const
{
	if (Size != a.Size || StartIndex != a.StartIndex)
		return false;
	for (int i = 0; i < Size; i++)
	{
		if (pVector[i] != a.pVector[i])
			return false;
	}
	return true;
}

template <class ValType> // ���������
bool TVector<ValType>::operator!=(const TVector &a) const
{
	if (Size != a.Size || StartIndex != a.StartIndex)
		return true;
	for (int i = 0; i < Size; i++)
	{
		if (pVector[i] != a.pVector[i])
			return true;
	}
	return false;
}

template <class ValType> // ������������
TVector<ValType>& TVector<ValType>::operator=(const TVector &a)
{
	if (&a != this)
	{
		if (Size != a.Size)
		{
			delete[] pVector;
			pVector = new ValType[a.Size];
			Size = a.Size;
		}
		StartIndex = a.StartIndex;
		for (int i = 0; i < Size; i++)
		{
			pVector[i] = a.pVector[i];
		}
	}
	return *this;
}

template <class ValType> // ��������� ������
TVector<ValType> TVector<ValType>::operator+(const ValType &val)
{
	TVector res(Size, StartIndex);
	for (int i = 0; i < Size; i++)
	{
		res.pVector[i] = pVector[i] + val;
	}
	return res;
}

template <class ValType> // ������� ������
TVector<ValType> TVector<ValType>::operator-(const ValType &val)
{
	TVector res(Size, StartIndex);
	for (int i = 0; i < Size; i++)
	{
		res.pVector[i] = pVector[i] - val;
	}
	return res;
}

template <class ValType> // �������� �� ������
TVector<ValType> TVector<ValType>::operator*(const ValType &val)
{
	TVector res(Size, StartIndex);
	for (int i = 0; i < Size; i++)
	{
		res.pVector[i] = pVector[i] * val;
	}
	return res;
}

template <class ValType> // �������� �� ������
TVector<ValType> TVector<ValType>::operator*(const int &val)
{
	TVector res(Size, StartIndex);
	for (int i = 0; i < Size; i++)
	{
		res.pVector[i] = pVector[i] * val;
	}
	return res;
}

template <class ValType> // ��������
TVector<ValType> TVector<ValType>::operator+(const TVector<ValType> &v)
{
	if (Size != v.Size)
		throw - 1;
	TVector res(Size, StartIndex);
	for (int i = 0; i < Size; i++)
	{
		res.pVector[i] = pVector[i] + v.pVector[i];
	}
	return res;
}

template <class ValType> // ���������
TVector<ValType> TVector<ValType>::operator-(const TVector<ValType> &v)
{
	if (Size != v.Size)
		throw - 1;
	TVector res(Size, StartIndex);
	for (int i = 0; i < Size; i++)
	{
		res.pVector[i] = pVector[i] - v.pVector[i];
	}
	return res;
}

template <class ValType> // ��������� ������������
ValType TVector<ValType>::operator*(const TVector<ValType> &v)
{
	if (Size != v.Size)
		throw - 1;
	ValType res = 0;
	for (int i = 0; i < Size; i++)
	{
		res += pVector[i] * v.pVector[i];
	}
	return res;
}
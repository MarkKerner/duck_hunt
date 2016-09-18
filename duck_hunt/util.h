#ifndef UTIL_H_MARK
#define UTIL_H_MARK


template <typename T>
void print_vector(vector<T> v)
{
	if (v.size() == 0)
		return;

	cerr << " print vector, with size " << v.size() << endl;
	cerr << v[0];
	for (int i = 1; i < v.size(); ++i)
	{
		cerr << " " << v[i];
	}
	cerr << endl;
}


#endif // !UTIL_H_MARK
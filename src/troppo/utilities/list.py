from numpy.core._multiarray_umath import ndarray, array


def is_list(x):
	return type(x) in [list, tuple] and len(x) > 0 or isinstance(x, ndarray) and x.size > 0


def is_list_else_empty(x):
	return type(x) in [list, tuple] or isinstance(x, ndarray)


def if_none_return_list(x):
	if x is None:
		return array([])
	else:
		return x
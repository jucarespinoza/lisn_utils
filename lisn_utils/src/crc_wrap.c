/*   crc_wrap.c
 *   Functions wich are too sloooow in native python.
 *
 *   Juan C. Espinoza, June 2011
 */

#define PROG_VERSION "1.0"

#include "Python.h"

#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif

extern unsigned int CRC32 (unsigned int, unsigned char *);

static PyObject *
lib_CRC32 (PyObject *self, PyObject *args)
{
  
  const char *ucBuffer;
  PyObject *input;
  Py_ssize_t ulCount;
  int ulCRC;

  if (! PyArg_ParseTuple (args, "O", &input))
	return NULL;

  if (PyObject_AsCharBuffer (input, &ucBuffer, &ulCount) != 0)
	return NULL;
  
  ulCRC = CRC32 (ulCount, (unsigned char *) ucBuffer);

  return Py_BuildValue ("i", ulCRC);

}


static PyMethodDef LibMethods[] = {
	{ (char *)"CRC32", lib_CRC32, METH_VARARGS, "Calculate CRC32 of 'char buffer' data"},
	{ NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC
init_crc(void)
{
  (void) Py_InitModule ("_crc", LibMethods);
}

#include <Python.h>

#include "extractor-core/include/extractor.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>


static PyObject*
extractor_extract_dna(PyObject*, PyObject* args)
{
    char* reference;
    char* sample;

    if (!PyArg_ParseTuple(args, "ss", &reference, &sample))
    {
        return NULL;
    } // if

    size_t const reference_length = strlen(reference);
    size_t const sample_length = strlen(sample);

    std::vector<mutalyzer::Variant> const variants = mutalyzer::extract_dna(reference, reference_length, sample, sample_length);

    PyObject* result = PyList_New(0);
    if (result == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_New");
        return NULL;
    } // if

    for (std::vector<mutalyzer::Variant>::const_iterator it = variants.begin(); it != variants.end(); ++it)
    {
        PyObject* item = Py_BuildValue("{s:i, s:i, s:i, s:i, s:i, s:i, s:i}",
                                       "reference_start", it->reference_start,
                                       "reference_end", it->reference_end,
                                       "sample_start", it->sample_start,
                                       "sample_end", it->sample_end,
                                       "type", it->type,
                                       "transposition_start", it->transposition_start,
                                       "transposition_end", it->transposition_end);
        if (item == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        if (PyList_Append(result, item) != 0)
        {
            Py_DECREF(item);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_Append");
            return NULL;
        } // if
    } // for
    return result;
} // extractor_extract


static PyMethodDef ExtractorMethods[] =
{
    {"extract_dna", extractor_extract_dna, METH_VARARGS, "Run the extract function."},

    {NULL, NULL, 0, NULL}  // sentinel
}; // ExtractorMethods


static struct PyModuleDef extractormodule =
{
    PyModuleDef_HEAD_INIT,
    "description-extractor",
    "HGVS Variant Description Extractor",
    -1,
    ExtractorMethods
}; // extractormodule


PyMODINIT_FUNC
PyInit_extractor(void)
{
    return PyModule_Create(&extractormodule);
} // PyInit_extractor


int
main(int, char* argv[])
{
    wchar_t* const program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL)
    {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        return EXIT_FAILURE;
    } // if

    PyImport_AppendInittab("description-extractor", PyInit_extractor);

    Py_SetProgramName(program);

    Py_Initialize();

    PyMem_RawFree(program);

    return EXIT_SUCCESS;
} // main

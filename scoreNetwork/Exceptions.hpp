#ifndef EXCEPTIONS_HPP_29045684568920568920456590268
#define EXCEPTIONS_HPP_29045684568920568920456590268

#include <exception>
#include <stdexcept>

class TypeException: public std::exception
{
    virtual const char* what() const throw()
    {
	return "Type exception";
    }
} ; //typeEx;

class GenericError: public std::runtime_error
{
public:
    GenericError(const std::string& msg = ""):std::runtime_error(msg) {}
} ;

#endif // EXCEPTIONS_HPP_


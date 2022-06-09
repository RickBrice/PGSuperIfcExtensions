#pragma once

// Use this exception class when we encounter IFC tags that can't be handled.
// Right now, it just gives a simple text message, but this class can be made
// more safisticated in the future if need be.
class CIfcAlignmentConverterException
{
public:
   CIfcAlignmentConverterException(LPCTSTR strMsg);
   ~CIfcAlignmentConverterException(void);

   LPCTSTR What() const;

private:
   std::_tstring m_strWhat;
};

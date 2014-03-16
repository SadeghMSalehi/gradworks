CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(muparser)

INCLUDE_DIRECTORIES(/opt/local/include)

SET (muparser_SRC 
muParserBase.cpp muParserBase.h muParserBytecode.cpp muParserBytecode.h muParserCallback.cpp muParserCallback.h muParser.cpp muParserDef.h muParserDLL.cpp muParserDLL.h muParserError.cpp muParserError.h muParserFixes.h muParser.h muParserInt.cpp muParserInt.h muParserStack.h muParserTemplateMagic.h muParserTest.cpp muParserTest.h muParserToken.h muParserTokenReader.cpp muParserTokenReader.h)

ADD_LIBRARY( muparser STATIC ${muparser_SRC})

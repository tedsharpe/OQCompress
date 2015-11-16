/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file BGZF.h
 * \author tsharpe
 * \date May 21, 2009
 *
 * \brief Utilities for writing BAM files.
 */
#ifndef LOOKUP_BGZF_H_
#define LOOKUP_BGZF_H_

#include <fstream>

class GZIPHeader
{
public:
    GZIPHeader()
    : mID1(31), mID2(139), mCompressionMethod(8), mFlags(4), mModTime(0),
      mExtraFlags(0), mOpSys(255), mXLen(6), mSI1('B'), mSI2('C'), mSLen(2)
    {}

    // compiler-supplied copying and destructor are OK

private:
    unsigned char mID1;
    unsigned char mID2;
    unsigned char mCompressionMethod;
    unsigned char mFlags;
    unsigned int mModTime;
    unsigned char mExtraFlags;
    unsigned char mOpSys;
    unsigned short mXLen;
    unsigned char mSI1;
    unsigned char mSI2;
    unsigned short mSLen;
protected:
    unsigned short mBlockSizeLessOne;
};

class BGZFBlock : public GZIPHeader
{
public:
    // compiler-supplied no-arg constructor, copying and destructor are OK

    // returns the amount of uncompressed data put into the block.
    // if len is too big, or the data is too incompressible, not all of the supplied
    // data will fit into a block.
    unsigned int compress( void* data, unsigned int len );

    unsigned int getBlockSize()
    { return mBlockSizeLessOne + 1U; }

private:
    bool tryCompress( void* data, unsigned int& len );

    // mBlockSizeLessOne is 16 bits, so 64K is the max block length and there are 26 bytes of header and footer
    unsigned char mDataBlock[64*1024UL-26UL];
    unsigned int mCRC32; // these last two members actually immediately follow the variable-length data block
    unsigned int mInputSize; // this struct is what a maximum-size block looks like

    static int const WINDOW_BITS = -15;
    static int const MEM_LEVEL = 8;
};

class BGZFStreambuf : public std::streambuf
{
public:
    BGZFStreambuf( std::streambuf* psb )
    : mpSB(psb)
    { setp(mBuf,mBuf+sizeof(mBuf)-1); }

    ~BGZFStreambuf()
    { sync(); }

private:
    BGZFStreambuf( BGZFStreambuf const& ); // undefined -- no copying
    BGZFStreambuf& operator=( BGZFStreambuf const& ); // undefined -- no copying

    int_type overflow( int_type ch );
    int sync();

    std::streambuf* mpSB;
    char mBuf[128UL*1024UL];
};

class BAMostream : public std::ostream
{
public:
    BAMostream( char const* bamFile );
    void close();

    std::filebuf mFilebuf;
    BGZFStreambuf mSB;
};

#endif /* LOOKUP_BGZF_H_ */

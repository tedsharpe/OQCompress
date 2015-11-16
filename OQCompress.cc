/*
 * OQCompress.cc
 *
 *  Created on: Nov 13, 2015
 *      Author: tsharpe
 */

#include "BGZF.h"
#include "gzstream.h"
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#define BAMERR(file,message)  \
     (std::cout << "\nBAM file " << file << message << '\n'), exit(1)

// image of the header for a BAM file
struct BAMAlignHead
{
    uint32_t mRemainingBlockSize;
    int32_t mRefID;
    int32_t mPos;
    uint8_t mNameLen;
    uint8_t mMapQ;
    uint16_t mBin;
    uint16_t mCigarLen;
    uint16_t mFlags;
    uint32_t mSeqLen;
    int32_t mMateRefID;
    int32_t mMatePos;
    int32_t mTLen;
};

// auxiliary tags signal the data type of the tag with these characters
// a return of 0 means "variable length"
// a return of -1 means "illegal data type specifier"
inline int getTagLength( char dataType )
{
    int tagLen;
    switch ( dataType )
    {
    case 'A': case 'c': case 'C': tagLen = 1; break;
    case 's': case 'S':           tagLen = 2; break;
    case 'i': case 'I': case 'f': tagLen = 4; break;
    case 'Z': case 'H': case 'B': tagLen = 0; break;
    default:                      tagLen = -1; break;
    }
    return tagLen;
}

// class to do quality score compression and decompression
class QualCompressor
{
public:
    QualCompressor() {}
    QualCompressor( QualCompressor const& )=delete;
    QualCompressor& operator=( QualCompressor const& )=delete;

    std::vector<uint8_t>& encode( std::vector<uint8_t> const& quals );
    std::vector<uint8_t>& decode( std::vector<uint8_t> const& packedQuals );

private:
    size_t packedSize() const
    { return std::accumulate(mBlocks.begin(),mBlocks.end(),0ul,
           []( size_t acc, Block const& blk ) { return acc+blk.size(); }); }

    void configureBlocks( std::vector<uint8_t> const& quals );

    static int nlz( uint32_t val )
    { return val ? __builtin_clz(val) : 32; }

    static int ceilLg2( uint32_t val )
    { return 32-nlz(val-1); }

    struct Block
    { Block( uint8_t nQs, uint8_t bits, uint8_t minQ )
      : mNQs(nQs), mBits(bits), mMinQ(minQ) {}
      unsigned size() const { return blockSize(mNQs,mBits); }
      static unsigned blockSize( unsigned nQs, unsigned nBits )
      { return (nQs*nBits+17+7)>>3; }
      uint8_t mNQs; uint8_t mBits; uint8_t mMinQ; };

    std::vector<Block> mBlocks;
    std::vector<unsigned> mCosts;
    std::vector<uint8_t> mBuffer;
};

void QualCompressor::configureBlocks( std::vector<uint8_t> const& quals )
{
    mBlocks.clear();
    mCosts.clear();
    mCosts.reserve(quals.size()+1);
    mCosts.push_back(0); // cost of an empty compressed qual vector
    auto beg = quals.begin();
    auto end = quals.end();
    auto itr = beg;
    while ( itr != end )
    {
        uint8_t const MAX_Q = 63;
        if ( *itr > MAX_Q )
        {   std::cout << "\nYour input reads are funny.  I found a quality score of "
                 << uint32_t(*itr) << ".\nThe maximum value that I allow is "
                 << uint32_t(MAX_Q) << ".\n" << std::endl;
            exit(1);    }

        auto iCost = mCosts.end();
        uint32_t minVal = uint32_t(*itr);
        uint32_t maxVal = minVal;
        uint32_t bits = 0;
        uint32_t prevCost = *--iCost;
        uint32_t nQs = 1;
        uint32_t bestCost = prevCost + Block::blockSize(nQs,bits);
        Block best(1,bits,minVal);
        auto itr2 = itr;
        ++itr;
        while ( itr2 != beg && nQs < 255 )
        {
            unsigned val = *--itr2;
            if ( val > maxVal )
                maxVal = val;
            if ( val < minVal )
                minVal = val;
            bits = ceilLg2(maxVal+1u-minVal);
            prevCost = *--iCost;
            unsigned curCost = prevCost + Block::blockSize(++nQs,bits);
            if ( curCost < bestCost )
            {
                bestCost = curCost;
                best = Block(nQs,bits,minVal);
            }
        }
        mCosts.push_back(bestCost);
        unsigned toRemove = best.mNQs - 1;
        if ( !toRemove )
            mBlocks.push_back(best);
        else
        {
            while ( toRemove > mBlocks.back().mNQs )
            {
                toRemove -= mBlocks.back().mNQs;
                mBlocks.pop_back();
            }
            if ( toRemove == mBlocks.back().mNQs )
                mBlocks.back() = best;
            else
            {
                mBlocks.back().mNQs -= toRemove;
                mBlocks.push_back(best);
            }
        }
    }
}

std::vector<uint8_t>& QualCompressor::encode( std::vector<uint8_t> const& quals )
{
    configureBlocks(quals);

    mBuffer.reserve(quals.size());
    mBuffer.clear();
    auto itr = quals.begin();
    for ( Block const& block : mBlocks )
    {
        uint64_t nQs = block.mNQs;
        uint64_t nBits = block.mBits;
        uint64_t minQ = block.mMinQ;
        mBuffer.push_back(nQs);
        uint64_t bits = nBits;
        bits |= minQ << 3;
        mBuffer.push_back(bits);
        bits >>= 8;
        if ( !nBits )
        {
            mBuffer.push_back(bits);
            itr += nQs;
        }
        else
        {
            uint64_t off = 1;
            while ( nQs-- )
            {
                uint64_t val = *itr - minQ;
                ++itr;
                bits |= val << off;
                if ( (off += nBits) >= 8 )
                {
                    mBuffer.push_back(bits);
                    off -= 8;
                    bits >>= 8;
                }
            }
            if ( off )
                mBuffer.push_back(bits);
        }
    }
    mBuffer.push_back(0);
    return mBuffer;
}

std::vector<uint8_t>& QualCompressor::decode( std::vector<uint8_t> const& packedQuals )
{
    mBuffer.reserve(4*packedQuals.size());
    mBuffer.clear();
    if ( packedQuals.empty() )
        return mBuffer;
    uint64_t addr = reinterpret_cast<uint64_t>(&packedQuals[0]);
    uint64_t* buf = reinterpret_cast<uint64_t*>(addr&~7);
    uint64_t bits = *buf++;
    addr = (addr & 7) << 3;
    uint64_t remain = 64 - addr;
    bits >>= addr;
    uint64_t nQs;
    while ( (nQs = bits&0xff) )
    {
        bits >>= 8;
        if ( !(remain -= 8) )
        {
            bits = *buf++;
            remain = 64;
        }
        uint64_t nBits = bits & 0x07;
        bits >>= 3;
        uint64_t minQ = bits & 0x3f;
        bits >>= 6;
        if ( remain < 9 )
        {
            bits = *buf++;
            minQ |= (bits & 1) << 5;
            bits >>= 1;
            remain += 64;
        }
        remain -= 9;
        if ( !nBits )
            while ( nQs-- )
                mBuffer.push_back(minQ);
        else
        {
            uint64_t mask = (1ul<<nBits)-1ul;
            while ( nQs-- )
            {
                uint64_t val = bits;
                uint64_t used = nBits;
                if ( remain < nBits )
                {
                    bits = *buf++;
                    val |= bits << remain;
                    used -= remain;
                    remain += 64;
                }
                remain -= nBits;
                bits >>= used;
                mBuffer.push_back(minQ + (val & mask));
            }
        }
        bits >>= remain & 7;
        remain &= ~7ul;
        if ( !remain )
        {
            bits = *buf++;
            remain = 64;
        }
    }
    return mBuffer;
}

int main( int argc, char** argv )
{
    if ( argc != 3 )
    {
        std::cout << "Usage: OQCompress in.bam out.bam" << std::endl;
        exit(1);
    }
    char const* inFile = argv[1];
    char const* outFile = argv[2];
    igzstream is(inFile);
    BAMostream os(outFile);

    // copy header
    std::vector<char> buffer;
    buffer.reserve(2048);
    uint32_t val;
    if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
        BAMERR(inFile," is empty");
    if ( val != 0x014d4142 )
        BAMERR(inFile," lacks a BAM header");
    if ( !os.write(reinterpret_cast<char const*>(&val),sizeof(val)) )
        BAMERR(outFile," is unwritable");
    if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
        BAMERR(inFile," header length is truncated");
    if ( !os.write(reinterpret_cast<char const*>(&val),sizeof(val)) )
        BAMERR(outFile," header length unwritable");
    buffer.resize(val);
    if ( !is.read(&buffer[0],val) )
        BAMERR(inFile," header is truncated");
    if ( !os.write(&buffer[0],val) )
        BAMERR(outFile," header unwritable");

    // copy reference dictionary
    uint32_t nRefs;
    if ( !is.read(reinterpret_cast<char*>(&nRefs),sizeof(nRefs)) )
        BAMERR(inFile," is truncated at ref desc count");
    if ( !os.write(reinterpret_cast<char const*>(&nRefs),sizeof(nRefs)) )
        BAMERR(outFile," ref desc count unwritable");
    while ( nRefs-- )
    {
        if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
            BAMERR(inFile," is truncated in ref desc len");
        if ( !os.write(reinterpret_cast<char*>(&val),sizeof(val)) )
            BAMERR(outFile," ref desc len unwritable");
        buffer.resize(val);
        if ( !is.read(&buffer[0],val) )
            BAMERR(inFile," ref desc name is truncated");
        if ( !os.write(&buffer[0],val) )
            BAMERR(outFile," ref desc name unwritable");
        if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
            BAMERR(inFile," is truncated in ref desc size");
        if ( !os.write(reinterpret_cast<char const*>(&val),sizeof(val)) )
            BAMERR(outFile," ref desc size unwritable");
    }

    BAMAlignHead aln;
    QualCompressor qc;
    size_t alnNo = 0;
    while ( is.peek() != std::istream::traits_type::eof() )
    {
        std::ostringstream oss;
        if ( !is.read(reinterpret_cast<char*>(&aln),sizeof(aln)) )
            BAMERR(inFile," is truncated in alignment header " << alnNo);
        buffer.resize(aln.mNameLen);
        if ( !is.read(&buffer[0],aln.mNameLen) )
            BAMERR(inFile," is truncated in read name " << alnNo);
        if ( !oss.write(&buffer[0],aln.mNameLen) )
            BAMERR(outFile," read name " << alnNo << " unwritable");
        val = aln.mCigarLen*sizeof(uint32_t);
        int32_t auxLen = aln.mRemainingBlockSize -
                                sizeof(aln) + sizeof(aln.mRemainingBlockSize) -
                                aln.mNameLen - val;
        buffer.resize(val);
        if ( !is.read(&buffer[0],val) )
            BAMERR(inFile," is truncated in cigar" << alnNo);
        if ( !oss.write(&buffer[0],val) )
            BAMERR(outFile," cigar " << alnNo << " unwritable");
        val = (aln.mSeqLen + 1)/2;
        auxLen -= val + aln.mSeqLen;
        if ( auxLen < 0 )
            BAMERR(inFile," invalid alignment block size" << alnNo);
        buffer.resize(val);
        if ( !is.read(&buffer[0],val) )
            BAMERR(inFile," sequence " << alnNo << " is truncated ");
        if ( !oss.write(&buffer[0],val) )
            BAMERR(outFile," sequence " << alnNo << " unwritable");
        buffer.resize(aln.mSeqLen);
        if ( !is.read(&buffer[0],aln.mSeqLen) )
            BAMERR(inFile," quals " << alnNo << " truncated");
        if ( !oss.write(&buffer[0],aln.mSeqLen) )
            BAMERR(outFile," quals " << alnNo << " unwritable");
        while ( auxLen > 0 )
        {
            char tag[3];
            if ( !is.read(tag,sizeof(tag)) )
                BAMERR(inFile," tag header truncated in alignment " << alnNo);
            if ( tag[0] == 'O' && tag[1] == 'Q' )
            {
                if ( tag[2] != 'Z' )
                    BAMERR(inFile," contains OQ tag with non-Z data type in alignment " << alnNo);
                buffer.resize(aln.mSeqLen);
                if ( !is.read(&buffer[0],aln.mSeqLen) )
                    BAMERR(inFile," is truncated in OQ tag data in alignment " << alnNo);
                char byte;
                if ( !is.get(byte) || byte )
                    BAMERR(inFile," contains OQ tag with the wrong length in alignment " << alnNo);

                for ( char& val : buffer )
                    val -= 33;

                std::vector<uint8_t> const& byteBuf = *reinterpret_cast<std::vector<uint8_t> const*>(&buffer);
                std::vector<uint8_t> const& packedQuals = qc.encode(byteBuf);
                if ( !oss.write("ZQBC",4) )
                    BAMERR(outFile," ZQ tag header unwritable in alignment " << alnNo);
                uint32_t size = packedQuals.size();
                if ( !oss.write(reinterpret_cast<char const*>(&size),sizeof(size)) )
                    BAMERR(outFile," ZQ tag length unwritable in alignment " << alnNo);
                if ( size && !oss.write(reinterpret_cast<char const*>(&packedQuals[0]),size) )
                    BAMERR(outFile," ZQ tag contents unwritable in alignment " << alnNo);

                auxLen -= aln.mSeqLen + 1 + 3;
                continue;
            }

            if ( tag[0] == 'Z' && tag[1] == 'Q' )
            {
                if ( tag[2] != 'B' )
                    BAMERR(inFile," contains a ZQ tag with non-B data type in alignment " << alnNo);
                char dataType;
                if ( !is.get(dataType) || dataType != 'C' )
                    BAMERR(inFile," contains a ZQ tag with non-C data type in alignment " << alnNo);
                uint32_t size;
                if ( !is.read(reinterpret_cast<char*>(&size),sizeof(size)) )
                    BAMERR(inFile," ZQ tag size truncated in alignment " << alnNo);
                buffer.resize(size);
                if ( !is.read(&buffer[0],size) )
                    BAMERR(inFile," ZQ tag data truncated in alignment " << alnNo);
                std::vector<uint8_t> const& byteBuf = *reinterpret_cast<std::vector<uint8_t> const*>(&buffer);
                std::vector<uint8_t>& quals = qc.decode(byteBuf);
                if ( quals.size() != aln.mSeqLen )
                    BAMERR(inFile," unpacked ZQ tag has wrong size in alignment " << alnNo);

                for ( uint8_t& val : quals )
                    val += 33;

                if ( !oss.write("OQZ",3) )
                    BAMERR(outFile," OQ tag header unwritable in alignment " << alnNo);
                if ( !oss.write(reinterpret_cast<char const*>(&quals[0]),quals.size()) )
                    BAMERR(outFile," OQ tag data unwritable in alignment " << alnNo);
                if ( !oss.put(0) )
                    BAMERR(outFile," OQ tag null unwritable in alignment " << alnNo);

                auxLen -= size + 4 + 4;
                continue;
            }

            if ( !oss.write(tag,sizeof(tag)) )
                BAMERR(outFile," tag header unwritable in alignment " << alnNo);
            auxLen -= 3;
            int tagLen = getTagLength(tag[2]);
            if ( tagLen == -1 )
                BAMERR(inFile," has bad data type in tag header in alignment " << alnNo);
            if ( tag[2] == 'B' )
            {
                char dataType;
                uint32_t arrLen;
                if ( !is.get(dataType) ||
                     !is.read(reinterpret_cast<char*>(&arrLen),sizeof(arrLen)) )
                    BAMERR(inFile," is truncated in B tag header in alignment " << alnNo);
                if ( !oss.put(dataType) ||
                     !oss.write(reinterpret_cast<char const*>(&arrLen),sizeof(arrLen)) )
                    BAMERR(outFile," B tag header unwritable in alignment " << alnNo);
                tagLen = getTagLength(dataType);
                if ( tagLen <= 0 )
                    BAMERR(inFile," has bad data type in B tag header in alignment " << alnNo);
                tagLen *= arrLen;
                auxLen -= 5;
            }

            if ( tagLen )
            {
                buffer.resize(tagLen);
                if ( !is.read(&buffer[0],tagLen) )
                    BAMERR(inFile," is truncated in tag data in alignment " << alnNo);
                if ( !oss.write(&buffer[0],tagLen) )
                    BAMERR(outFile," tag data in alignment " << alnNo << " unwritable");
                auxLen -= tagLen;
            }
            else // must be H or Z tag type
            {
                char byte;
                do
                {
                    if ( !is.get(byte) )
                        BAMERR(inFile," is truncated in null-delimited tag data for alignment " << alnNo);
                    if ( !oss.put(byte) )
                        BAMERR(outFile," null-delimited tag data in alignment " << alnNo << " unwritable");
                    auxLen -= 1;
                } while ( byte );
            }
        }

        if ( auxLen < 0 )
            BAMERR(inFile," has bogus alignment block len for alignment " << alnNo);
        std::string alnBytes = oss.str();
        aln.mRemainingBlockSize = alnBytes.size() + sizeof(aln) - sizeof(aln.mRemainingBlockSize);
        if ( !os.write(reinterpret_cast<char const*>(&aln),sizeof(aln)) )
            BAMERR(outFile," alignment header in alignment " << alnNo);
        if ( !os.write(&alnBytes[0],alnBytes.size()) )
            BAMERR(outFile," alignment data in alignment " << alnNo);
    }
}

// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "rpcserver.h"
#include "main.h"
#include "kernel.h"
#include "checkpoints.h"

using namespace json_spirit;
using namespace std;

extern void TxToJSON(const CTransaction& tx, const uint256 hashBlock, json_spirit::Object& entry);
extern enum Checkpoints::CPMode CheckpointsMode;

double GetDifficulty(const CBlockIndex* blockindex, int algo, bool next)
{
    // Floating point number that is a multiple of the minimum difficulty,
    // minimum difficulty = 1.0.
    if (blockindex == NULL)
    {
        if (pindexBest == NULL)
            return 1.0;
	else if (algo<0)
	  blockindex = GetLastBlockIndex(pindexBest, false);
	else
	  blockindex = pindexBest;
    }

    if (algo<0) {
      if (blockindex->IsProofOfStake()) {
	algo = 1;
      }
      else {
	algo = 0;
      }
    }
    unsigned int nBits = 0;
    bool fProofOfStake = false;
    if (algo == 1) {
      fProofOfStake = true;
    }
    if (next) {
      nBits = GetNextTargetRequired(blockindex,fProofOfStake);
    }
    else {
      bool pos_tip = blockindex->IsProofOfStake();
      if (pos_tip != fProofOfStake) {
	blockindex = GetLastBlockIndex(blockindex,fProofOfStake);
      }
      nBits = blockindex->nBits;
    }
    
    int nShift = (nBits >> 24) & 0xff;

    double dDiff =
        (double)0x0000ffff / (double)(nBits & 0x00ffffff);

    while (nShift < 29)
    {
        dDiff *= 256.0;
        nShift++;
    }
    while (nShift > 29)
    {
        dDiff /= 256.0;
        nShift--;
    }

    return dDiff;
}

double GetPeakHashrate (const CBlockIndex* blockindex, int algo) {
  const int nBlocksDay = 720;
  if (blockindex == NULL) {
      if (pindexBest == NULL)
	return 0.;
      else
	blockindex = pindexBest;
  }
  bool pos = false;
  if (algo == 1) pos = true;
  bool posTip = blockindex->IsProofOfStake();
  int algoTip = posTip ? 1 : 0;
  if (algoTip != algo) {
    blockindex = GetPrevBlockIndex(blockindex,algo);
  }
  if (!blockindex) return 0.;
  do {
    if ((blockindex->nHeight % 1440) == 0) {
      double hashes_peak = 0.;
      const CBlockIndex * pprev_algo = GetPrevBlockIndex(blockindex,algo);
      for (int i=0; i<365; i++) {
	if (!pprev_algo) break;
	int time_f = pprev_algo->GetMedianTimePast();
	CBigNum hashes_bn = CBigNum(pprev_algo->GetBlockTrust());
	int time_i = 0;
	
	for (int j=0; j<nBlocksDay-1; j++) {
	 
	  pprev_algo = GetPrevBlockIndex(pprev_algo,-1);

	  if (pprev_algo) {
	    time_i = pprev_algo->GetMedianTimePast();
	  }
	  else {
	    hashes_bn = CBigNum(0);
	    break;
	  }
	  //LogPrintf("j=%d add block work of block %lu\n",j,pprev_algo->nHeight);
	  hashes_bn += CBigNum(pprev_algo->GetBlockTrust());
	}
	const CBlockIndex * pprev_algo_time = GetPrevBlockIndex(pprev_algo,-1);
	if (pprev_algo_time) {
	  time_i = pprev_algo_time->GetMedianTimePast();
	}
	else {
	  const CBlockIndex * blockindex_time = pprev_algo;
	  while (blockindex_time && blockindex_time->nHeight>=1) {
	    blockindex_time = blockindex_time->pprev;
	  }
	  if (blockindex_time) {
	    time_i = blockindex_time->GetBlockTime();
	  }
	}
	pprev_algo = pprev_algo_time;
	
	if (time_f>time_i) {
	  time_f -= time_i;
	}
	else {
	  return std::numeric_limits<double>::max();
	}
	//LogPrintf("hashes = %f, time = %f\n",(double)hashes_bn.getulong(),(double)time_f);
	double hashes = (hashes_bn/time_f).getuint256().getdouble();
	//LogPrintf("hashes per sec = %f\n",hashes);
	if (hashes>hashes_peak) hashes_peak = hashes;
      }
      return hashes_peak;
      break;
    }
    blockindex = blockindex->pprev;
  } while (blockindex);
  return 0.;
}

double GetCurrentHashrate (const CBlockIndex* blockindex, int algo) {
  const int nBlocksDay = 720;
  if (blockindex == NULL)
    {
      if (pindexBest == NULL)
	return 0.;
      else
	blockindex = pindexBest;
    }
  
  bool pos = false;
  if (algo == 1) pos = true;
  bool pos_tip = blockindex->IsProofOfStake();
  int algo_tip = pos_tip ? 1 : 0;
  if (algo_tip != algo) {
    blockindex = GetPrevBlockIndex(blockindex,pos);
  }
  if (!blockindex) {
    return 0.;
  }
  do {
    if ((blockindex->nHeight % 1440) == 0) {
      const CBlockIndex * pcur_algo = GetPrevBlockIndex(blockindex,algo);
      if (!pcur_algo) return 0.;
      int time_f = pcur_algo->GetMedianTimePast();
      CBigNum hashes_bn = CBigNum(pcur_algo->GetBlockTrust());
      int time_i = 0;
      const CBlockIndex * pprev_algo = pcur_algo;
      for (int j=0; j<nBlocksDay-1; j++) {
	pprev_algo = GetPrevBlockIndex(pprev_algo,-1);
	if (pprev_algo) {
	  time_i = pprev_algo->GetMedianTimePast();
	}
	else {
	  return 0.;
	}
	hashes_bn += CBigNum(pprev_algo->GetBlockTrust());
      }
      const CBlockIndex * pprev_algo_time = GetPrevBlockIndex(pprev_algo,-1);
      if (pprev_algo_time) {
	time_i = pprev_algo_time->GetMedianTimePast();
      }
      else {
	const CBlockIndex * blockindex_time = pprev_algo;
	while (blockindex_time && blockindex_time->nHeight>=1) {
	  blockindex_time = blockindex_time->pprev;
	}
	if (blockindex_time) time_i = blockindex_time->GetBlockTime();
      }

      if (time_f>time_i) {
	time_f -= time_i;
      }
      else {
	return std::numeric_limits<double>::max();
      }
      return (hashes_bn/time_f).getuint256().getdouble();
    }
    blockindex = blockindex->pprev;
  } while (blockindex);
  return 0.;
}

double GetAverageBlockSpacing (const CBlockIndex * blockindex, const int algo, const int averagingInterval) {
  
  if (averagingInterval <= 1) return 0.;

  if (blockindex == NULL) {
    if (pindexBest == NULL)
      return 0.;
    else
      blockindex = pindexBest;
  }
  
  const CBlockIndex *BlockReading = blockindex;
  int64_t CountBlocks = 0;
  int64_t nActualTimespan = 0;
  int64_t LastBlockTime = 0;

  bool pos = false;
  if (algo == 1) pos = true;
  
  for (unsigned int i = 1; BlockReading && BlockReading->nHeight > 0; i++) {
    if (CountBlocks >= averagingInterval) { break; }
    bool block_pos = BlockReading->IsProofOfStake();
    if (algo >=0 && block_pos != pos) {
      BlockReading = BlockReading->pprev;
      continue;
    }
    CountBlocks++;
    if(LastBlockTime > 0){
      nActualTimespan = LastBlockTime - BlockReading->GetMedianTimePast();
    }
    else {
      LastBlockTime = BlockReading->GetMedianTimePast();
    }

    BlockReading = BlockReading->pprev;
    
  }
  return ((double)nActualTimespan)/((double)averagingInterval)/60.;
}

double GetPoWMHashPS()
{
    if (pindexBest->nHeight >= LAST_POW_BLOCK)
        return 0;

    int nPoWInterval = 72;
    int64_t nTargetSpacingWorkMin = 30, nTargetSpacingWork = 30;

    CBlockIndex* pindex = pindexGenesisBlock;
    CBlockIndex* pindexPrevWork = pindexGenesisBlock;

    while (pindex)
    {
        if (pindex->IsProofOfWork())
        {
            int64_t nActualSpacingWork = pindex->GetBlockTime() - pindexPrevWork->GetBlockTime();
            nTargetSpacingWork = ((nPoWInterval - 1) * nTargetSpacingWork + nActualSpacingWork + nActualSpacingWork) / (nPoWInterval + 1);
            nTargetSpacingWork = max(nTargetSpacingWork, nTargetSpacingWorkMin);
            pindexPrevWork = pindex;
        }

        pindex = pindex->pnext;
    }

    return GetDifficulty() * 4294.967296 / nTargetSpacingWork;
}

double GetNetworkHashPS(int lookup, int height, int algo) {
  if (lookup == -1) return GetPoWMHashPS();
    CBlockIndex *pb = pindexBest;

    if (height >= 0 && height <pindexBest->nHeight) {
      while (pb != 0 && pb->nHeight > height) {
	pb = pb->pprev;
      }
    }

    if (pb == NULL || !pb->nHeight)
        return 0;

    // If lookup is -1, then use blocks since last difficulty change.
    if (lookup <= 0)
        lookup = pb->nHeight % 2016 + 1;

    // If lookup is larger than chain, then set it to chain length.
    if (lookup > pb->nHeight)
        lookup = pb->nHeight;

    const CBlockIndex *pb0 = pb;
    bool pos = false;
    if (algo == 1) pos = true;
    bool pos_tip = pb0->IsProofOfStake();
    if (pos_tip != pos) {
      pb0 = GetPrevBlockIndex(pb0,algo);
    }
    if (!pb0) return 0.;
    int64_t minTime = pb0->GetBlockTime();
    int64_t maxTime = minTime;
    CBigNum hashes_bn = CBigNum(pb0->GetBlockTrust());
    for (int i = 0; i < lookup; i++) {
        pb0 = pb0->pprev;
	if (!pb0) break;
	if (pb0->IsProofOfStake()!=pos) {
	  lookup++;
	  continue;
	}
        int64_t time = pb0->GetBlockTime();
        minTime = std::min(time, minTime);
        maxTime = std::max(time, maxTime);
	hashes_bn += CBigNum(pb0->GetBlockTrust());
    }

    // In case there's a situation where minTime == maxTime, we don't want a divide by zero exception.
    if (minTime == maxTime)
        return 0;

    //uint256 workDiff = pb->nChainWork - pb0->nChainWork;
    int64_t timeDiff = maxTime - minTime;

    return ((hashes_bn.getuint256().getdouble()) / ((double)timeDiff));
}

double GetPoSKernelPS()
{
    int nPoSInterval = 72;
    double dStakeKernelsTriedAvg = 0;
    int nStakesHandled = 0, nStakesTime = 0;

    CBlockIndex* pindex = pindexBest;;
    CBlockIndex* pindexPrevStake = NULL;

    while (pindex && nStakesHandled < nPoSInterval)
    {
        if (pindex->IsProofOfStake())
        {
            if (pindexPrevStake)
            {
                dStakeKernelsTriedAvg += GetDifficulty(pindexPrevStake) * 4294967296.0;
                nStakesTime += pindexPrevStake->nTime - pindex->nTime;
                nStakesHandled++;
            }
            pindexPrevStake = pindex;
        }

        pindex = pindex->pprev;
    }

    double result = 0;

    if (nStakesTime){
        result = dStakeKernelsTriedAvg / nStakesTime;
        result *= STAKE_TIMESTAMP_MASK + 1;
    }

    return result;
}

Object blockToJSON(const CBlock& block, const CBlockIndex* blockindex, bool fPrintTransactionDetail)
{
    Object result;
    result.push_back(Pair("hash", block.GetHash().GetHex()));
    int confirmations = -1;
    // Only report confirmations if the block is on the main chain
    if (blockindex->IsInMainChain())
        confirmations = nBestHeight - blockindex->nHeight + 1;
    result.push_back(Pair("confirmations", confirmations));
    result.push_back(Pair("size", (int)::GetSerializeSize(block, SER_NETWORK, PROTOCOL_VERSION)));
    result.push_back(Pair("height", blockindex->nHeight));
    result.push_back(Pair("version", block.nVersion));
    result.push_back(Pair("merkleroot", block.hashMerkleRoot.GetHex()));
    result.push_back(Pair("mint", ValueFromAmount(blockindex->nMint)));
    result.push_back(Pair("time", (int64_t)block.GetBlockTime()));
    result.push_back(Pair("mediantime", (int64_t)(blockindex->GetMedianTimePast())));
    result.push_back(Pair("nonce", (uint64_t)block.nNonce));
    result.push_back(Pair("bits", strprintf("%08x", block.nBits)));
    result.push_back(Pair("difficulty", GetDifficulty(blockindex)));
    result.push_back(Pair("blocktrust", leftTrim(blockindex->GetBlockTrust().GetHex(), '0')));
    result.push_back(Pair("chaintrust", leftTrim(blockindex->nChainTrust.GetHex(), '0')));
    if (blockindex->pprev)
        result.push_back(Pair("previousblockhash", blockindex->pprev->GetBlockHash().GetHex()));
    if (blockindex->pnext)
        result.push_back(Pair("nextblockhash", blockindex->pnext->GetBlockHash().GetHex()));

    result.push_back(Pair("flags", strprintf("%s%s", blockindex->IsProofOfStake()? "proof-of-stake" : "proof-of-work", blockindex->GeneratedStakeModifier()? " stake-modifier": "")));
    result.push_back(Pair("proofhash", blockindex->hashProof.GetHex()));
    result.push_back(Pair("entropybit", (int)blockindex->GetStakeEntropyBit()));
    result.push_back(Pair("modifier", strprintf("%016x", blockindex->nStakeModifier)));
    Array txinfo;
    BOOST_FOREACH (const CTransaction& tx, block.vtx)
    {
        if (fPrintTransactionDetail)
        {
            Object entry;

            entry.push_back(Pair("txid", tx.GetHash().GetHex()));
            TxToJSON(tx, 0, entry);

            txinfo.push_back(entry);
        }
        else
            txinfo.push_back(tx.GetHash().GetHex());
    }

    result.push_back(Pair("tx", txinfo));

    if (block.IsProofOfStake())
        result.push_back(Pair("signature", HexStr(block.vchBlockSig.begin(), block.vchBlockSig.end())));

    return result;
}

Value getbestblockhash(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getbestblockhash\n"
            "Returns the hash of the best block in the longest block chain.");

    return hashBestChain.GetHex();
}

Value getblockcount(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getblockcount\n"
            "Returns the number of blocks in the longest block chain.");

    return nBestHeight;
}

Value getdifficulty (const Array& params, bool fHelp) {
  if (fHelp || params.size() > 3)
    throw runtime_error(
			"getdifficulty ( algo height next )\n"
			"Returns an object containing difficulty info.\n"
				    "\nArguments:\n"
	    "1. \"algo\"     (numeric, optional) The algo, (miningalgo) by default\n"
	    "2. \"height\"     (numeric, optional) The height to look at, tip by default\n"
	    "3. \"next\"     (boolean, optional) Whether to get next difficulty required, false by default\n"
	    "\nResult:\n"
	    "{\n"
	    " \"difficulty\": xxxxx           (numeric)\n"
	    "}\n"
			);

  int algo = -1;
  CBlockIndex * blockindex = NULL;
  bool next = false;

  bool v2 = false;
  if (params.size()>0) {
    v2 = true;
    algo = params[0].get_int();
    if (params.size()>1) {
      int height = params[1].get_int();
      blockindex = pindexBest;
      if (height >= 0) {
	while (blockindex && blockindex->nHeight > height) {
	  blockindex = blockindex->pprev;
	}
      }
      if (params.size()>2) {
	next = params[2].get_bool();
      }
    }
  }

  Object obj;
  if (v2) {
    obj.push_back(Pair("difficulty",(double)GetDifficulty(blockindex,algo,next)));
  }
  else {
    obj.push_back(Pair("proof-of-work",        GetDifficulty()));
    obj.push_back(Pair("proof-of-stake",       GetDifficulty(GetLastBlockIndex(pindexBest, true))));
  }
  return obj;
}

Value getrawmempool(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getrawmempool\n"
            "Returns all transaction ids in memory pool.");

    vector<uint256> vtxid;
    mempool.queryHashes(vtxid);

    Array a;
    BOOST_FOREACH(const uint256& hash, vtxid)
        a.push_back(hash.ToString());

    return a;
}

Value getblockhash(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "getblockhash <index>\n"
            "Returns hash of block in best-block-chain at <index>.");

    int nHeight = params[0].get_int();
    if (nHeight < 0 || nHeight > nBestHeight)
        throw runtime_error("Block number out of range.");

    CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
    return pblockindex->phashBlock->GetHex();
}

Value getblock(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 2)
        throw runtime_error(
            "getblock <hash> [txinfo]\n"
            "txinfo optional to print more detailed tx info\n"
            "Returns details of a block with given block-hash.");

    std::string strHash = params[0].get_str();
    uint256 hash(strHash);

    if (mapBlockIndex.count(hash) == 0)
        throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Block not found");

    CBlock block;
    CBlockIndex* pblockindex = mapBlockIndex[hash];
    block.ReadFromDisk(pblockindex, true);

    return blockToJSON(block, pblockindex, params.size() > 1 ? params[1].get_bool() : false);
}

Value getblockbynumber(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 2)
        throw runtime_error(
            "getblockbynumber <number> [txinfo]\n"
            "txinfo optional to print more detailed tx info\n"
            "Returns details of a block with given block-number.");

    int nHeight = params[0].get_int();
    if (nHeight < 0 || nHeight > nBestHeight)
        throw runtime_error("Block number out of range.");

    CBlock block;
    CBlockIndex* pblockindex = mapBlockIndex[hashBestChain];
    while (pblockindex->nHeight > nHeight)
        pblockindex = pblockindex->pprev;

    uint256 hash = *pblockindex->phashBlock;

    pblockindex = mapBlockIndex[hash];
    block.ReadFromDisk(pblockindex, true);

    return blockToJSON(block, pblockindex, params.size() > 1 ? params[1].get_bool() : false);
}

// ppcoin: get information of sync-checkpoint
Value getcheckpoint(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "getcheckpoint\n"
            "Show info of synchronized checkpoint.\n");

    Object result;
    CBlockIndex* pindexCheckpoint;

    result.push_back(Pair("synccheckpoint", Checkpoints::hashSyncCheckpoint.ToString()));
    pindexCheckpoint = mapBlockIndex[Checkpoints::hashSyncCheckpoint];
    result.push_back(Pair("height", pindexCheckpoint->nHeight));
    result.push_back(Pair("timestamp", DateTimeStrFormat(pindexCheckpoint->GetBlockTime())));

    // Check that the block satisfies synchronized checkpoint
    if (CheckpointsMode == Checkpoints::STRICT)
        result.push_back(Pair("policy", "strict"));

    if (CheckpointsMode == Checkpoints::ADVISORY)
        result.push_back(Pair("policy", "advisory"));

    if (CheckpointsMode == Checkpoints::PERMISSIVE)
        result.push_back(Pair("policy", "permissive"));

    if (mapArgs.count("-checkpointkey"))
        result.push_back(Pair("checkpointmaster", true));

    return result;
}

Value chaindynamics(const Array& params, bool fHelp)
{
    if (fHelp || params.size() > 1)
        throw runtime_error(
            "chaindynamics (height)\n"
            "Returns an object containing various state info.\n"
            "}\n"
	    "\nResult:\n"
	    "{\n"
	    " \"difficulty <algo>\": xxxxx           (numeric),\n"
	    " \"peak hashrate <algo>\": xxxxx           (numeric),\n"
	    " \"current hashrate <algo>\": xxxxx           (numeric),\n"
	    " \"average block spacing <algo>\": xxxxx           (numeric)\n"
	    "}\n"
        );

    CBlockIndex * pindex = 0;
    if (params.size()>0) {
      int height = params[0].get_int();
      pindex = pindexBest;
      while (pindex && pindex->nHeight > height)
        pindex = pindex->pprev;
    }    
    
    Object obj;

    obj.push_back(Pair("difficulty SHA256D",    (double)GetDifficulty(pindex,0,true)));
    obj.push_back(Pair("difficulty PoS",    (double)GetDifficulty(pindex,1,true)));
    obj.push_back(Pair("peak hashrate SHA256D",    (double)GetPeakHashrate(pindex,0)));
    obj.push_back(Pair("peak hashrate PoS",    (double)GetPeakHashrate(pindex,1)));
    obj.push_back(Pair("current hashrate SHA256D",    (double)GetCurrentHashrate(pindex,0)));    
    obj.push_back(Pair("current hashrate PoS",    (double)GetCurrentHashrate(pindex,1)));
    obj.push_back(Pair("average block spacing",    (double)GetAverageBlockSpacing(pindex,-1)));
    obj.push_back(Pair("average block spacing SHA256D",    (double)GetAverageBlockSpacing(pindex,0)));    
    obj.push_back(Pair("average block spacing PoS",    (double)GetAverageBlockSpacing(pindex,1)));

    return obj;
}

Value getblockspacing(const Array& params, bool fHelp)
{
    if (fHelp)
        throw runtime_error(
            "getblockspacing (algo interval height )\n"
            "Returns an object containing blockspacing info.\n"
	    "\nArguments:\n"
	    "1. \"algo\"     (numeric, optional) The algo (-1) by default\n"
            "2. \"interval\"     (numeric, optional) The interval in number of blocks, 25 by default\n"
	    "3. \"height\"     (numeric, optional) The height for the endpoint of the interval, tip by default\n"	    
	    "\nResult:\n"
	    "{\n"
	    "  \"average block spacing\": xxxxx           (numeric)\n"
	    "}\n"
			    );

    int algo = -1;
    int interval = 25;
    CBlockIndex * blockindex = NULL;
    
    if (params.size()>0) {
      algo = params[0].get_int();
      if (params.size()>1) {
	interval = params[1].get_int();
	if (params.size()>2) {
	  int height = params[2].get_int();
	  blockindex = pindexBest;
	  while (blockindex && blockindex->nHeight > height)
	    blockindex = blockindex->pprev;
	}
      }
    }
    
    Object obj;
    obj.push_back(Pair("average block spacing",    (double)GetAverageBlockSpacing(blockindex,algo,interval)));

    return obj;
}

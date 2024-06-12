// Original Code: Copyright (c) 2011-2014 The Bitcoin Core Developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#define BOOST_TEST_MODULE eMark Test Suite



#include "main.h"
#include "txdb.h"
#include "ui_interface.h"
#include "util.h"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>


extern bool fPrintToConsole;
extern void noui_connect();

struct TestingSetup {
  //CCoinsViewDB *pcoinsdbview;
    boost::filesystem::path pathTemp;
    boost::thread_group threadGroup;

    TestingSetup() {
      fPrintToDebugLog = false; // don't want to write to debug.log file
        noui_connect();

        pathTemp = boost::filesystem::temp_directory_path() / strprintf("test_eMark_%lu_%i", (unsigned long)GetTime(), (int)(GetRand(100000)));
        boost::filesystem::create_directories(pathTemp);
        mapArgs["-datadir"] = pathTemp.string();
	SelectParams(CChainParams::MAIN);
        //pblocktree = new CBlockTreeDB(1 << 20, true);
        //pcoinsdbview = new CCoinsViewDB(1 << 23, true);
        //pcoinsTip = new CCoinsViewCache(*pcoinsdbview);
        //InitBlockIndex();
        //nScriptCheckThreads = 3;
        //for (int i=0; i < nScriptCheckThreads-1; i++)
	//  threadGroup.create_thread(&ThreadScriptCheck);
        RegisterNodeSignals(GetNodeSignals());
    }
    ~TestingSetup()
    {
        threadGroup.interrupt_all();
        threadGroup.join_all();
        UnregisterNodeSignals(GetNodeSignals());
        //delete pcoinsTip;
        //delete pcoinsdbview;
        //delete pblocktree;
        boost::filesystem::remove_all(pathTemp);
    }
};

BOOST_GLOBAL_FIXTURE(TestingSetup);

void Shutdown(void* parg)
{
  exit(0);
}

/*void StartShutdown()
{
  exit(0);
  }

bool ShutdownRequested()
{
  return false;
  }
*/

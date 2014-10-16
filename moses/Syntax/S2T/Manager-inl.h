#pragma once

#include "moses/DecodeGraph.h"
#include "moses/StaticData.h"
#include "moses/Syntax/BoundedPriorityContainer.h"
#include "moses/Syntax/CubeQueue.h"
#include "moses/Syntax/PHyperedge.h"
#include "moses/Syntax/RuleTable.h"
#include "moses/Syntax/RuleTableFF.h"
#include "moses/Syntax/SHyperedgeBundle.h"
#include "moses/Syntax/SVertex.h"
#include "moses/Syntax/SVertexRecombinationOrderer.h"
#include "moses/Syntax/SymbolEqualityPred.h"
#include "moses/Syntax/SymbolHasher.h"

#include "OovHandler.h"
#include "PChart.h"
#include "RuleTrie.h"
#include "SChart.h"

namespace Moses
{
namespace Syntax
{
namespace S2T
{

template<typename Parser>
Manager<Parser>::Manager(const InputType &source)
    : m_source(source)
{
}

template<typename Parser>
void Manager<Parser>::InitializePChart(PChart &pchart)
{
  // Create cells.
  pchart.cells.resize(m_source.GetSize());
  for (std::size_t i = 0; i < m_source.GetSize(); ++i) {
    pchart.cells[i].resize(m_source.GetSize());
  }
  // Insert PVertex objects for source words
  for (std::size_t i = 0; i < m_source.GetSize(); ++i) {
    PVertex v(WordsRange(i,i), m_source.GetWord(i));
    PChart::Cell::TMap::value_type x(v.symbol, v);
    pchart.cells[i][i].terminalVertices.insert(x);
  }
}

template<typename Parser>
void Manager<Parser>::InitializeSChart(const PChart &pchart, SChart &schart)
{
  // Create cells.
  schart.cells.resize(m_source.GetSize());
  for (std::size_t i = 0; i < m_source.GetSize(); ++i) {
    schart.cells[i].resize(m_source.GetSize());
  }
  // Insert SVertex objects for source words
  for (std::size_t i = 0; i < m_source.GetSize(); ++i) {
    const Word &terminal = m_source.GetWord(i);
    boost::shared_ptr<SVertex> v(new SVertex());
    v->best = 0;
    const PChart::Cell::TMap &pmap = pchart.cells[i][i].terminalVertices;
    PChart::Cell::TMap::const_iterator p = pmap.find(terminal);
    assert(p != pmap.end());
    v->pvertex = &(p->second);
    SVertexBeam beam(1, v);
    SChart::Cell::TMap::value_type x(terminal, beam);
    schart.cells[i][i].terminalBeams.insert(x);
  }
}

template<typename Parser>
void Manager<Parser>::InitializeParsers(PChart &pchart,
                                        std::size_t ruleLimit)
{
  const std::vector<RuleTableFF*> &ffs = RuleTableFF::Instances();

  const std::vector<DecodeGraph*> &graphs =
    StaticData::Instance().GetDecodeGraphs();

  assert(ffs.size() == graphs.size());

  for (std::size_t i = 0; i < ffs.size(); ++i) {
    RuleTableFF *ff = ffs[i];
    std::size_t maxChartSpan = graphs[i]->GetMaxChartSpan();
    // This may change in the future, but currently we assume that every
    // RuleTableFF is associated with a static, file-based rule table of
    // some sort and that the table should have been loaded into a RuleTable
    // by this point.
    const RuleTable *table = ff->GetTable();
    assert(table);
    RuleTable *nonConstTable = const_cast<RuleTable*>(table);
    boost::shared_ptr<Parser> parser;
    typename Parser::RuleTrie *trie =
      dynamic_cast<typename Parser::RuleTrie*>(nonConstTable);
    if (!trie) {
      // FIXME
      assert(false);
    }
    parser.reset(new Parser(pchart, *trie, maxChartSpan));
    m_parsers.push_back(parser);
  }

  // Check for OOVs and synthesize an additional rule trie + parser if
  // necessary.
  std::set<Word> oovs;
  std::size_t maxOovWidth = 0;
  FindOovs(pchart, oovs, maxOovWidth);
  if (!oovs.empty()) {
    // FIXME Add a hidden RuleTableFF for unknown words(?)
    OovHandler<typename Parser::RuleTrie> oovHandler(*ffs[0]);
    m_oovRuleTrie = oovHandler.SynthesizeRuleTrie(oovs.begin(), oovs.end());
    // Create a parser for the OOV rule trie.
    boost::shared_ptr<Parser> parser(
        new Parser(pchart, *m_oovRuleTrie, maxOovWidth));
    m_parsers.push_back(parser);
  }
}

// Find the set of OOVs for this input.  This function assumes that the
// PChart argument has already been initialized from the input.
template<typename Parser>
void Manager<Parser>::FindOovs(const PChart &pchart, std::set<Word> &oovs,
                               std::size_t maxOovWidth)
{
  // Get the set of RuleTries.
  std::vector<const RuleTrie *> tries;
  const std::vector<RuleTableFF*> &ffs = RuleTableFF::Instances();
  for (std::size_t i = 0; i < ffs.size(); ++i) {
    const RuleTableFF *ff = ffs[i];
    if (ff->GetTable()) {
      const RuleTrie *trie = dynamic_cast<const RuleTrie*>(ff->GetTable());
      assert(trie);  // FIXME
      tries.push_back(trie);
    }
  }

  // For every sink vertex in pchart (except for <s> and </s>), check whether
  // the word has a preterminal rule in any of the rule tables.  If not then
  // add it to the OOV set.
  oovs.clear();
  maxOovWidth = 0;
  // Assume <s> and </s> have been added at sentence boundaries, so skip
  // cells starting at position 0 and ending at the last position.
  for (std::size_t i = 1; i < pchart.cells.size()-1; ++i) {
    for (std::size_t j = i; j < pchart.cells[i].size()-1; ++j) {
      std::size_t width = j-i+1;
      const PChart::Cell::TMap &map = pchart.cells[i][j].terminalVertices;
      for (PChart::Cell::TMap::const_iterator p = map.begin();
           p != map.end(); ++p) {
        const Word &word = p->first;
        assert(!word.IsNonTerminal());
        bool found = false;
        for (std::vector<const RuleTrie *>::const_iterator q = tries.begin();
             q != tries.end(); ++q) {
          const RuleTrie *trie = *q;
          if (trie->HasPreterminalRule(word)) {
            found = true;
            break;
          }
        }
        if (!found) {
          oovs.insert(word);
          maxOovWidth = std::max(maxOovWidth, width);
        }
      }
    }
  }
}

template<typename Parser>
void Manager<Parser>::Decode()
{
  const StaticData &staticData = StaticData::Instance();

  // Get various pruning-related constants.
  const std::size_t popLimit = staticData.GetCubePruningPopLimit();
  const std::size_t ruleLimit = staticData.GetRuleLimit();
  const std::size_t beamLimit = staticData.GetMaxHypoStackSize();

  // Initialise PChart and SChart
  InitializePChart(m_pchart);
  InitializeSChart(m_pchart, m_schart);

  // Initialize the parsers.
  InitializeParsers(m_pchart, ruleLimit);

  // Create a callback to process the PHyperedges produced by the parsers.
  typename Parser::CallbackType callback(m_schart, ruleLimit);

  // Visit each cell of PChart in right-to-left bottom-up order.
  std::size_t size = m_source.GetSize();
  for (int start = size-1; start >= 0; --start) {
    for (std::size_t width = 1; width <= size-start; ++width) {
      std::size_t end = start + width - 1;

      PChart::Cell &pcell = m_pchart.cells[start][end];
      SChart::Cell &scell = m_schart.cells[start][end];

      WordsRange range(start, end);

      // Call the parsers to generate PHyperedges for this span and convert
      // each one to a SHyperedgeBundle (via the callback).  The callback
      // prunes the SHyperedgeBundles and keeps the best ones (up to ruleLimit).
      callback.InitForRange(range);
      for (typename std::vector<boost::shared_ptr<Parser> >::iterator
           p = m_parsers.begin(); p != m_parsers.end(); ++p) {
        (*p)->EnumerateHyperedges(range, callback);
      }

      // Retrieve the (pruned) set of SHyperedgeBundles from the callback.
      const BoundedPriorityContainer<SHyperedgeBundle> &bundles =
          callback.GetContainer();

      // Use cube pruning to extract SHyperedges from SHyperedgeBundles.
      // Collect the SHyperedges into buffers, one for each category.
      CubeQueue cubeQueue(bundles.Begin(), bundles.End());
      std::size_t count = 0;
      typedef boost::unordered_map<Word, std::vector<SHyperedge*>,
                                   SymbolHasher, SymbolEqualityPred > BufferMap;
      BufferMap buffers;
      while (count < popLimit && !cubeQueue.IsEmpty()) {
        SHyperedge *hyperedge = cubeQueue.Pop();
        // BEGIN{HACK}
        // The way things currently work, the LHS of each hyperedge is not
        // determined until just before the point of its creation, when a
        // target phrase is selected from the list of possible phrases (which
        // happens during cube pruning).  The cube pruning code doesn't (and
        // shouldn't) know about the contents of PChart and so creation of
        // the PVertex is deferred until this point.
        const Word &lhs = hyperedge->translation->GetTargetLHS();
        hyperedge->head->pvertex =
            pcell.nonTerminalVertices.Insert(lhs, PVertex(range, lhs));
        // END{HACK}
        buffers[lhs].push_back(hyperedge);
        ++count;
      }

      // Recombine SVertices and sort into beams.
      for (BufferMap::const_iterator p = buffers.begin(); p != buffers.end();
           ++p) {
        const Word &category = p->first;
        const std::vector<SHyperedge*> &buffer = p->second;
        // FIXME
        SVertexBeam *beam = scell.nonTerminalBeams.Find(category);
        if (!beam) {
          beam = scell.nonTerminalBeams.Insert(category, SVertexBeam());
        }
        RecombineAndSort(buffer, *beam);
      }

      // Prune beams.
      if (beamLimit > 0) {
        for (SChart::Cell::NMap::Iterator p = scell.nonTerminalBeams.Begin();
             p != scell.nonTerminalBeams.End(); ++p) {
          SVertexBeam &beam = p->second;
          if (beam.size() > beamLimit) {
            beam.resize(beamLimit);
          }
        }
      }

      // Prune the PChart cell for this span by removing vertices for
      // categories that don't occur in the SChart.
// Note: see HACK above.  Pruning the chart isn't currently necessary.
//      PrunePChart(scell, pcell);
    }
  }
}

template<typename Parser>
const SHyperedge *Manager<Parser>::GetBestSHyperedge() const
{
  const SChart::Cell::NMap &beams =
      m_schart.cells[0][m_source.GetSize()-1].nonTerminalBeams;
  if (beams.Size() == 0) {
    return 0;
  }
  assert(beams.Size() == 1);
  const std::vector<boost::shared_ptr<SVertex> > &beam = beams.Begin()->second;
  return beam[0]->best;
}

template<typename Parser>
void Manager<Parser>::ExtractKBest(
    std::size_t k,
    std::vector<boost::shared_ptr<KBestExtractor::Derivation> > &kBestList,
    bool onlyDistinct) const
{
  kBestList.clear();
  if (k == 0 || m_source.GetSize() == 0) {
    return;
  }

  // Get the top-level SVertex beam.
  const SChart::Cell::NMap &beams =
      m_schart.cells[0][m_source.GetSize()-1].nonTerminalBeams;
  if (beams.Size() == 0) {
    return;
  }
  assert(beams.Size() == 1);
  const std::vector<boost::shared_ptr<SVertex> > &beam = beams.Begin()->second;

  KBestExtractor extractor;

  if (!onlyDistinct) {
    // Return the k-best list as is, including duplicate translations.
    extractor.Extract(beam, k, kBestList);
    return;
  }

  // Determine how many derivations to extract.  If the k-best list is
  // restricted to distinct translations then this limit should be bigger
  // than k.  The k-best factor determines how much bigger the limit should be,
  // with 0 being 'unlimited.'  This actually sets a large-ish limit in case
  // too many translations are identical.
  const StaticData &staticData = StaticData::Instance();
  const std::size_t nBestFactor = staticData.GetNBestFactor();
  std::size_t numDerivations = (nBestFactor == 0) ? k*1000 : k*nBestFactor;

  // Extract the derivations.
  KBestExtractor::KBestVec bigList;
  bigList.reserve(numDerivations);
  extractor.Extract(beam, numDerivations, bigList);

  // Copy derivations into kBestList, skipping ones with repeated translations.
  std::set<Phrase> distinct;
  for (KBestExtractor::KBestVec::const_iterator p = bigList.begin();
       kBestList.size() < k && p != bigList.end(); ++p) {
    boost::shared_ptr<KBestExtractor::Derivation> derivation = *p;
    Phrase translation = KBestExtractor::GetOutputPhrase(*derivation);
    if (distinct.insert(translation).second) {
      kBestList.push_back(derivation);
    }
  }
}

template<typename Parser>
void Manager<Parser>::PrunePChart(const SChart::Cell &scell,
                                  PChart::Cell &pcell)
{
/* FIXME
  PChart::Cell::VertexMap::iterator p = pcell.vertices.begin();
  while (p != pcell.vertices.end()) {
    const Word &category = p->first;
    if (scell.beams.find(category) == scell.beams.end()) {
      PChart::Cell::VertexMap::iterator q = p++;
      pcell.vertices.erase(q);
    } else {
      ++p;
    }
  }
*/
}

template<typename Parser>
void Manager<Parser>::RecombineAndSort(const std::vector<SHyperedge*> &buffer,
                                       SVertexBeam &beam)
{
  // Step 1: Create a map containing a single instance of each distinct vertex
  // (where distinctness is defined by the state value).  The hyperedges'
  // head pointers are updated to point to the vertex instances in the map and
  // any 'duplicate' vertices are deleted.
// TODO Set?
  typedef std::map<SVertex *, SVertex *, SVertexRecombinationOrderer> Map;
  Map map;
  for (std::vector<SHyperedge*>::const_iterator p = buffer.begin();
       p != buffer.end(); ++p) {
    SHyperedge *h = *p;
    SVertex *v = h->head;
    assert(v->best == h);
    assert(v->recombined.empty());
    std::pair<Map::iterator, bool> result = map.insert(Map::value_type(v, v));
    if (result.second) {
      continue;  // v's recombination value hasn't been seen before.
    }
    // v is a duplicate (according to the recombination rules).
    // Compare the score of h against the score of the best incoming hyperedge
    // for the stored vertex.
    SVertex *storedVertex = result.first->second;
    if (h->score > storedVertex->best->score) {
      // h's score is better.
      storedVertex->recombined.push_back(storedVertex->best);
      storedVertex->best = h;
    } else {
      storedVertex->recombined.push_back(h);
    }
    h->head->best = 0;
    delete h->head;
    h->head = storedVertex;
  }

  // Step 2: Copy the vertices from the map to the beam.
  beam.clear();
  beam.reserve(map.size());
  for (Map::const_iterator p = map.begin(); p != map.end(); ++p) {
    beam.push_back(boost::shared_ptr<SVertex>(p->first));
  }

  // Step 3: Sort the vertices in the beam.
  std::sort(beam.begin(), beam.end(), SVertexBeamContentOrderer());
}

}  // S2T
}  // Syntax
}  // Moses
/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_integration.hh  --  ONE PLACE WHERE THE NEW XML BLOCKS ARE READ AND ACTED ON
 * ================================================================================================
 *
 *  The five modules added in this round -- SBML reading, model fetching, geometry generation,
 *  surrogate training and CSV diagnostics -- are each self-contained and each testable on their
 *  own. This header is the seam between them and complab.cpp.
 *
 *  IT EXISTS SO THAT complab.cpp CHANGES AS LITTLE AS POSSIBLE. That file is 150 kB and cannot be
 *  compiled outside Palabos, so every line added to it is a line that cannot be tested until the
 *  first cluster build. Everything that CAN live outside it, does: all the XML parsing, all the
 *  validation, all the messages. complab.cpp gains five short calls.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE FIVE CALL SITES, in the order complab.cpp reaches them
 *
 *    1. after initialize_complab()          integ::readConfig(cfg)
 *    2. before the geometry is read         integ::provideGeometry(cfg, ...)
 *    3. while loading metabolic models      integ::resolveModel(cfg, name, ...)
 *    4. after the models are loaded         integ::prepareSurrogate(cfg, ...)
 *    5. at each save interval               integ::recordDiagnostics(cfg, ...)
 *
 *  Each is a no-op unless the matching XML block is present, so an existing input file behaves
 *  exactly as it did before.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT IS NOT TESTED
 *
 *  Everything in the five modules is tested standalone. This header is NOT: it needs Palabos to
 *  compile, so the first time it is built is on the cluster. It is written to be as close to
 *  trivial as possible for that reason -- parsing, string building, and calls into code that has
 *  been exercised. There is no arithmetic here worth getting wrong.
 * ================================================================================================
 */
#ifndef COMPLAB3D_INTEGRATION_HH
#define COMPLAB3D_INTEGRATION_HH

#include <string>
#include <vector>

#include "complab3d_sbml.hh"
#include "complab3d_fetch.hh"
#include "complab3d_geometry.hh"
#include "complab3d_surrogate.hh"
#include "complab3d_diagnostics.hh"

namespace integ {

/* ------------------------------------------------------------------------------------------------
 *  Everything the new blocks configure, in one struct.
 * ------------------------------------------------------------------------------------------------ */
struct Config {
    /* <model_source> and friends */
    std::string modelSource;          // "bigg:iJO1366", or empty
    std::string modelCache;           // default "input"
    std::string modelBundle;          // default "models"
    bool allowDownload;

    /* <domain><generate> / <import_raw> */
    bool generateGeometry;
    complab_geom::GenOptions gen;
    std::string importRaw;
    std::string rawDtype;
    double rawThreshold;
    bool rawInvert;
    std::string writeGeometryTo;      // where to put the generated .dat, for the record

    /* <surrogate> */
    bool srgEnabled;
    std::string srgWeights;
    bool srgTrainIfMissing;
    /* Which microbe this network is for. -1, the default, means every microbe whose
     * <reaction_type> is a surrogate type -- the common case of one organism, where making the
     * user write an index they cannot get wrong is just an opportunity to get it wrong. */
    int srgMicrobe;
    complab_srg::TrainOptions srgTrain;

    /* <diagnostics> */
    bool diagEnabled;
    std::string diagCsv;
    long diagInterval;
    double diagTol;
    std::vector<std::string> diagConserve;

    Config()
      : modelCache("input"), modelBundle("models"), allowDownload(false),
        generateGeometry(false), rawDtype("uint8"), rawThreshold(128.0), rawInvert(false),
        srgEnabled(false), srgTrainIfMissing(false), srgMicrobe(-1),
        diagEnabled(false), diagInterval(0), diagTol(1e-6) {}
};

/* ------------------------------------------------------------------------------------------------
 *  Parsing helpers.
 *
 *  Written against Palabos's XMLreader through a tiny shim so that this header does not need
 *  Palabos to be READ by a human, and so the pattern "absent means keep the default" is written
 *  once rather than forty times. Every optional tag behaves the same way: absent is fine, present
 *  but malformed stops the run with the tag named.
 * ------------------------------------------------------------------------------------------------ */
template <class Reader, class T>
inline bool getOpt(Reader &doc, const char *a, const char *b, const char *c, T &out)
{
    try {
        if (c) doc[a][b][c].read(out);
        else if (b) doc[a][b].read(out);
        else doc[a].read(out);
        return true;
    } catch (...) { return false; }
}

inline bool truthy(std::string s)
{
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] = (char) (s[i] - 'A' + 'a');
    return s == "true" || s == "yes" || s == "1" || s == "on";
}

/* Split "acetate 1e-3 10 log" into an input specification. */
inline bool parseInputSpec(const std::string &line, complab_srg::InputSpec &s, std::string &err)
{
    std::vector<std::string> tok;
    std::string cur;
    for (size_t i = 0; i <= line.size(); ++i) {
        if (i == line.size() || line[i] == ' ' || line[i] == '\t') {
            if (!cur.empty()) { tok.push_back(cur); cur.clear(); }
        } else cur += line[i];
    }
    if (tok.size() < 3) {
        err = "'" + line + "' should read: <substrate name> <low> <high> [log]";
        return false;
    }
    s.name = tok[0];
    s.lo = std::atof(tok[1].c_str());
    s.hi = std::atof(tok[2].c_str());
    s.logScale = (tok.size() > 3 && (tok[3] == "log" || tok[3] == "logarithmic"));
    if (s.hi <= s.lo) { err = "'" + line + "': the high end must exceed the low end"; return false; }
    if (s.logScale && s.lo <= 0.0) {
        err = "'" + line + "': a logarithmic sweep needs a strictly positive low end";
        return false;
    }
    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 1: read the new blocks. Returns "" on success, or a message to print before stopping.
 * ------------------------------------------------------------------------------------------------ */
template <class Reader>
inline std::string readConfig(Reader &doc, Config &cfg)
{
    const char *P = "parameters";
    std::string tmp;

    /* ---- model source ---- */
    if (getOpt(doc, P, "model_source", (const char *) 0, tmp)) cfg.modelSource = tmp;
    if (getOpt(doc, P, "model_cache", (const char *) 0, tmp)) cfg.modelCache = tmp;
    if (getOpt(doc, P, "model_bundle", (const char *) 0, tmp)) cfg.modelBundle = tmp;
    if (getOpt(doc, P, "allow_download", (const char *) 0, tmp)) cfg.allowDownload = truthy(tmp);

    /* ---- geometry ---- */
    if (getOpt(doc, "parameters", "LB_numerics", "domain", tmp)) { /* presence only */ }
    try {
        std::string kind;
        doc[P]["LB_numerics"]["domain"]["generate"].read(kind);
        cfg.generateGeometry = true;
        cfg.gen.kind = kind;
        double d;
        int i;
        if (getOpt(doc, P, "LB_numerics", "domain", tmp)) { }
        try { doc[P]["LB_numerics"]["domain"]["porosity"].read(d); cfg.gen.porosity = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["grain_radius"].read(d); cfg.gen.radius = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["grain_count"].read(i); cfg.gen.count = i; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["aperture"].read(d); cfg.gen.aperture = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["roughness"].read(d); cfg.gen.roughness = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["layers"].read(i); cfg.gen.layers = i; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["seed"].read(i); cfg.gen.seed = (unsigned long long) i; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["walls"].read(tmp); cfg.gen.walls = tmp; } catch (...) {}
    } catch (...) { cfg.generateGeometry = false; }

    try { doc[P]["LB_numerics"]["domain"]["import_raw"].read(cfg.importRaw); } catch (...) {}
    try { doc[P]["LB_numerics"]["domain"]["raw_dtype"].read(cfg.rawDtype); } catch (...) {}
    try { doc[P]["LB_numerics"]["domain"]["threshold"].read(cfg.rawThreshold); } catch (...) {}
    try { std::string s; doc[P]["LB_numerics"]["domain"]["invert"].read(s); cfg.rawInvert = truthy(s); } catch (...) {}
    try { doc[P]["LB_numerics"]["domain"]["write_geometry"].read(cfg.writeGeometryTo); } catch (...) {}

    if (cfg.generateGeometry && !cfg.importRaw.empty())
        return "  [GEOM] both <generate> and <import_raw> are set. Pick one.\n";

    /* ---- surrogate ---- */
    try {
        std::string en;
        doc[P]["surrogate"]["enabled"].read(en);
        cfg.srgEnabled = truthy(en);
    } catch (...) { cfg.srgEnabled = false; }

    if (cfg.srgEnabled) {
        try { doc[P]["surrogate"]["weights_file"].read(cfg.srgWeights); } catch (...) {}
        try { std::string s; doc[P]["surrogate"]["train_if_missing"].read(s); cfg.srgTrainIfMissing = truthy(s); } catch (...) {}
        try { int m; doc[P]["surrogate"]["microbe"].read(m); cfg.srgMicrobe = m; } catch (...) {}
        if (cfg.srgWeights.empty())
            return "  [SRG] <surrogate> is enabled but <weights_file> is not set.\n";

        int i;
        try { doc[P]["surrogate"]["train"]["samples"].read(i); cfg.srgTrain.samples = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["restarts"].read(i); cfg.srgTrain.restarts = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["epochs"].read(i); cfg.srgTrain.epochs = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["verify_points"].read(i); cfg.srgTrain.verifyPoints = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["seed"].read(i); cfg.srgTrain.seed = (unsigned long long) i; } catch (...) {}
        try { std::string s; doc[P]["surrogate"]["train"]["log_output"].read(s); cfg.srgTrain.logOutput = truthy(s); } catch (...) {}
        try {
            std::vector<int> L;
            doc[P]["surrogate"]["train"]["layers"].read(L);
            if (!L.empty()) cfg.srgTrain.hidden = L;
        } catch (...) {}

        cfg.srgTrain.inputs.clear();
        for (int k = 0; k < 8; ++k) {
            char tag[32];
            std::snprintf(tag, sizeof(tag), "input%d", k);
            std::string line;
            try { doc[P]["surrogate"]["train"][tag].read(line); } catch (...) { break; }
            complab_srg::InputSpec s;
            std::string err;
            if (!parseInputSpec(line, s, err))
                return "  [SRG] <" + std::string(tag) + ">: " + err + "\n";
            cfg.srgTrain.inputs.push_back(s);
        }
        if (cfg.srgTrainIfMissing && cfg.srgTrain.inputs.empty())
            return "  [SRG] <train_if_missing> is true but no <inputN> ranges are given, so there\n"
                   "  [SRG] is nothing to sweep. Add one line per input: <input0>name lo hi log</input0>\n";
        if (cfg.srgTrain.inputs.size() > 3)
            return "  [SRG] more than 3 training inputs. The sample count needed grows faster than\n"
                   "  [SRG] the fit improves; use 1 to 3, or fit offline with surrogate_training/.\n";
    }

    /* ---- diagnostics ---- */
    try {
        std::string en;
        doc[P]["diagnostics"]["enabled"].read(en);
        cfg.diagEnabled = truthy(en);
    } catch (...) { cfg.diagEnabled = false; }

    if (cfg.diagEnabled) {
        cfg.diagCsv = "summary.csv";
        try { doc[P]["diagnostics"]["summary_csv"].read(cfg.diagCsv); } catch (...) {}
        try { int i; doc[P]["diagnostics"]["interval"].read(i); cfg.diagInterval = i; } catch (...) {}
        try { doc[P]["diagnostics"]["tolerance"].read(cfg.diagTol); } catch (...) {}
        for (int k = 0; k < 32; ++k) {
            std::string e;
            if (k == 0) { try { doc[P]["diagnostics"]["conserve"].read(e); } catch (...) { break; } }
            else {
                char tag[32];
                std::snprintf(tag, sizeof(tag), "conserve%d", k);
                try { doc[P]["diagnostics"][tag].read(e); } catch (...) { break; }
            }
            if (!e.empty()) cfg.diagConserve.push_back(e);
        }
    }
    return "";
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 2: geometry.
 *
 *  Returns "" if the run should read <filename> as before. Otherwise generates or imports, writes
 *  a .dat so the run is reproducible from its own output, and returns the path to read.
 * ------------------------------------------------------------------------------------------------ */
inline std::string provideGeometry(const Config &cfg, int nx, int ny, int nz,
                                   const std::string &inputDir, std::string &log, bool &fatal)
{
    fatal = false;
    if (!cfg.generateGeometry && cfg.importRaw.empty()) return "";

    complab_geom::Grid g;
    std::string err;

    if (cfg.generateGeometry) {
        g = complab_geom::generate(nx, ny, nz, cfg.gen, &err);
        if (g.size() == 0) { log += "  [GEOM] " + err + "\n"; fatal = true; return ""; }
        log += "  [GEOM] generated a '" + cfg.gen.kind + "' pore space\n";
    } else {
        const std::string src = inputDir + cfg.importRaw;
        if (!complab_geom::importRaw(g, src, nx, ny, nz, cfg.rawDtype, cfg.rawThreshold,
                                     cfg.rawInvert, cfg.gen.walls, &err)) {
            log += "  [GEOM] " + err + "\n";
            fatal = true;
            return "";
        }
        log += "  [GEOM] imported " + src + "\n";
    }

    const complab_geom::Inspection R = complab_geom::inspect(g);
    log += complab_geom::inspectionText(g, R);
    if (!R.percolates) fatal = true;      // a sealed domain is not a well-posed flow problem

    const std::string out = inputDir + (cfg.writeGeometryTo.empty()
                                        ? std::string("generated_geometry.dat")
                                        : cfg.writeGeometryTo);
    if (complab_geom::writeDat(g, out, &err))
        log += "  [GEOM] written to " + out + " so this run can be reproduced exactly\n";
    else
        log += "  [GEOM] note: " + err + "\n";
    return out;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 3: find a model file, and say plainly whether it is the one expected.
 * ------------------------------------------------------------------------------------------------ */
inline std::string resolveModel(const Config &cfg, const std::string &fallbackName,
                                bool isMaster, std::string &log, bool &fatal)
{
    fatal = false;
    if (cfg.modelSource.empty())
        return cfg.modelCache + "/" + fallbackName + ".xml";

    complab_fetch::Result r = complab_fetch::ensureModel(cfg.modelSource, cfg.modelCache,
                                                         cfg.allowDownload, isMaster,
                                                         cfg.modelBundle);
    log += r.message;
    if (!r.ok) fatal = true;
    return r.path;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 4: load the surrogate, training it first if it is missing and training is allowed.
 * ------------------------------------------------------------------------------------------------ */
inline bool prepareSurrogate(const Config &cfg, complab_srg::Network &net,
                             complab_srg::FbaFn fba, void *ctx,
                             bool isMaster, std::string &log)
{
    if (!cfg.srgEnabled) return true;

    std::string err;
    if (complab_srg::load(net, cfg.srgWeights, &err)) {
        log += "  [SRG] loaded " + cfg.srgWeights + "\n";

        /* [FIX] Does this file belong to THIS model?
         *
         * Nothing used to check. Change <model_source> from bigg:e_coli_core to
         * bigg:iJO1366, forget to delete the old .srg, and the run quietly keeps using
         * the network fitted to the first organism -- every growth rate wrong, every
         * number plausible, and no line in the log to suggest it.
         *
         * The file records what it was trained on. Compare, and say so if they differ.
         * A warning rather than a stop: a network trained offline against the same model
         * under a different name is legitimate, and refusing to run would be worse than
         * saying what was noticed. */
        if (!cfg.modelSource.empty() && !net.provenance.empty()
            && net.provenance.find(cfg.modelSource) == std::string::npos) {
            log += "  [SRG] WARNING: this weights file says it was '" + net.provenance + "',\n"
                   "  [SRG] but <model_source> is '" + cfg.modelSource + "'. If the metabolic\n"
                   "  [SRG] model changed, this network was fitted to the OTHER one and every\n"
                   "  [SRG] growth rate it returns is wrong. Delete " + cfg.srgWeights + "\n"
                   "  [SRG] to refit, or ignore this if you know the two are the same model.\n";
        }
        char b[256];
        for (int i = 0; i < net.nIn; ++i) {
            std::snprintf(b, sizeof(b), "  [SRG]   input %d valid over %.6g .. %.6g\n",
                          i, net.trainMin[(size_t) i], net.trainMax[(size_t) i]);
            log += b;
        }
        log += "  [SRG] outside that box the network extrapolates and its answer means nothing.\n";
        return true;
    }

    if (!cfg.srgTrainIfMissing) {
        log += "  [SRG] " + err + "\n"
               "  [SRG] Set <train_if_missing>true</train_if_missing> to fit one now, or produce\n"
               "  [SRG] the file offline with surrogate_training/.\n";
        return false;
    }
    if (!fba) {
        log += "  [SRG] training needs a flux balance solver, and this executable was built\n"
               "  [SRG] without one. Rebuild with -DENABLE_GLPK=ON, or supply the weights file.\n";
        return false;
    }

    if (!isMaster) {
        /* Ranks other than 0 must not train: they would each fit a different network. The caller
         * barriers and then loads the file rank 0 wrote. */
        return true;
    }

    log += "  [SRG] no weights file; training one now. This happens once.\n";
    complab_srg::TrainReport rep;
    if (!complab_srg::train(net, cfg.srgTrain, fba, ctx, rep, &err)) {
        log += "  [SRG] training failed: " + err + "\n";
        return false;
    }
    net.provenance = "trained in-run from " + cfg.modelSource;
    log += complab_srg::reportText(net, rep);
    if (complab_srg::save(net, cfg.srgWeights, &err))
        log += "  [SRG] written to " + cfg.srgWeights + "; later runs will load it directly.\n";
    else
        log += "  [SRG] note: could not write " + cfg.srgWeights + " (" + err + ")\n";
    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 5: diagnostics. The caller has already reduced the totals across ranks.
 * ------------------------------------------------------------------------------------------------ */
inline void setupDiagnostics(const Config &cfg, complab_diag::Diagnostics &D, bool isMaster,
                             const std::string &outputDir = std::string())
{
    if (!cfg.diagEnabled) return;

    /* Put the CSV where the rest of the run's output goes. Without this it lands in whatever
     * directory the job was launched from, which on a cluster is neither where the user looks
     * nor, for an array job, unique. An absolute path in <summary_csv> is honoured as given. */
    std::string path = cfg.diagCsv;
    if (!outputDir.empty() && !path.empty() && path[0] != '/')
        path = outputDir + path;

    D.configure(true, path, cfg.diagTol, isMaster);
    for (size_t i = 0; i < cfg.diagConserve.size(); ++i) D.addConserve(cfg.diagConserve[i]);
}

}  // namespace integ

#endif  // COMPLAB3D_INTEGRATION_HH

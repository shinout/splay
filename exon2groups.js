var LineStream = require('linestream');
var color      = require('termcolor').define();

function main() {
  var filename = process.argv[2];
    console.assert(require("path").existsSync(filename));
  
  var stream = new LineStream(filename, {trim : true});

  var splices = {};
  var groups  = {};

  stream.on('data', function(line) {
    var vals = line.split('\t');
    if (vals.length < 11) return;

    var splicenames = vals[9].split(",").map(function(v) {
      return v.split("[")[0];
    });

    splicenames.forEach(function(spname) {
      if (!splices[spname]) {
        splices[spname] = new SpliceData(spname);
      }
    });

    var len = splicenames.length;

    for (var i=0; i < len-1; i++) {
      for (var j=0; j < len; j++) {
        var sp1 = splices[splicenames[i]];
        var sp2 = splices[splicenames[j]];
        sp1.relate(sp2);
      }
    }
  });

  stream.on('end', function() {
    Object.keys(splices).forEach(function(spname, k) {
      var sp = splices[spname];
      sp.visit(k + 1, function(sp, gid) {
        if (! groups[gid]) groups[gid] = [];
        groups[gid].push(sp.name);
      });
    });

    Object.keys(groups).forEach(function(k) {
      console.log(groups[k].join("\t"));
    });
  });
}

function SpliceData(name) {
  this.name = name;
  this.rels = {};
  this.gid   = null;
}

SpliceData.prototype.relate = function visit(spdata, noSetToTheOther) {
  this.rels[spdata.name] = spdata;
  if (!noSetToTheOther) spdata.relate(this, true);
};

SpliceData.prototype.visit = function(gid, fn) {
  if (this.gid != null) return;
  this.gid = gid;
  fn(this, gid);
  Object.keys(this.rels).forEach(function(spdname) {
    var spd = this.rels[spdname];
    spd.visit(gid, fn);
  }, this);
};


module.exports = main;
if (process.argv[1] === __filename) { main() }

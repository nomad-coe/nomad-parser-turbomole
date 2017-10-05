package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object TurbomoleParser extends SimpleExternalParserGenerator(
  name = "TurbomoleParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("TurbomoleParser")) ::
      ("parserId" -> jn.JString("TurbomoleParser" + lab.TurbomoleVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.TurbomoleVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s*(?<progr>[a-zA-z0-9_]+)\s*(?:\([^()]+\))\s*:\s*TURBOMOLE\s*(?<version>.*)
\s*Copyright \(C\) [0-9]+ TURBOMOLE GmbH, Karlsruhe
""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/turbomole/parser/parser-turbomole/TurbomoleParser.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-turbomole/TurbomoleESCFparser.py",
    "parser-turbomole/TurbomoleParser.py",
    "parser-turbomole/TurbomoleControlInParser.py",
    "parser-turbomole/TurbomoleCommon.py",
    "parser-turbomole/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/turbomole.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-turbomole" -> "parsers/turbomole/parser/parser-turbomole",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)

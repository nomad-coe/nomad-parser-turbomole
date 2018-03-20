/*
  Copyright 2016-2018 Arvid Conrad Ihrig, Aliaksei Mazheika
                      Fritz-Haber-Institut der Max-Planck-Gesellschaft

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
 */
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
    "parser-turbomole/EmbeddingParser.py",
    "parser-turbomole/GradientParser.py",
    "parser-turbomole/MethodParser.py",
    "parser-turbomole/OrbitalParser.py",
    "parser-turbomole/SystemParser.py",
    "parser-turbomole/AOFORCEparser.py",
    "parser-turbomole/CCSDF12parser.py",
    "parser-turbomole/DSCFparser.py",
    "parser-turbomole/ESCFparser.py",
    "parser-turbomole/GRADparser.py",
    "parser-turbomole/RICC2parser.py",
    "parser-turbomole/RIDFTparser.py",
    "parser-turbomole/RIRPAparser.py",
    "parser-turbomole/STATPTparser.py",
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

/*
 * Copyright 2016-2018 Sebastian Alarcon Villaseca, Sebastián Alarcón Villaseca, Fawzi Mohamed, Micael Oliveira, Ankit Kariryaa, Danio Brambila
 * 
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object OrcaParser extends SimpleExternalParserGenerator(
  name = "OrcaParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("OrcaParser")) ::
      ("parserId" -> jn.JString("OrcaParser" + lab.OrcaVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.OrcaVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*
\s+\* O   R   C   A \*\s*
\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*
\s*
\s*--- An Ab Initio, DFT and Semiempirical electronic structure package ---\s*
""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/orca/parser/parser-orca/orca_parser.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-orca/orca_parser.py",
    "parser-orca/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta.nomadmetainfo.json",
    "nomad_meta_info/orca.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-orca" -> "parsers/orca/parser/parser-orca",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping(),
  metaInfoEnv = Some(lab.meta.KnownMetaInfoEnvs.orca)
)

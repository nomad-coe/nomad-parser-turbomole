/*
   Copyright 2016-2017 The NOMAD Developers Group

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

import java.nio.file.Paths

import org.scalatest.{ Matchers, WordSpec }

class TurbomoleParserSpec extends WordSpec with Matchers {
  "the Turbomole Parser" when {
    val testFiles = Paths.get("test", "examples").toFile.listFiles().filter(_.toString.endsWith(".out"))
    val basePath = Paths.get("parsers", "turbomole")
    testFiles.foreach { file =>
      val filePath = basePath.resolve(file.toPath)
      "parsing '" + file.toString + "'" should {
        "return parsing results as json-events" in {
          val result = ParserRun.parse(TurbomoleParser, filePath.toString, "json-events")
          result should be(ParseResult.ParseSuccess)
        }
        "return parsing results as json" in {
          val result = ParserRun.parse(TurbomoleParser, filePath.toString, "json")
          result should be(ParseResult.ParseSuccess)
        }
      }
    }
  }
}

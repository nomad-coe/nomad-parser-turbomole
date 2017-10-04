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

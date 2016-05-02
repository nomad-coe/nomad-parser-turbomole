package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object TurbomoleParserSpec extends Specification {
  "TurbomoleParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(TurbomoleParser, "parsers/turbomole/test/examples", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(TurbomoleParser, "parsers/turbomole/test/examples", "json") must_== ParseResult.ParseSuccess
    }
  }
}

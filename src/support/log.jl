function log()

    println("\n", "="^99)

    println("$(StackTraces.stacktrace()[2].func)\n")

    println(now())

    println("="^99, "\n")

end
